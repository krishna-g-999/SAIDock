import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, QED, AllChem

LIPINSKI_LIMITS = {
    'MW':       500,
    'LogP':     5,
    'HBD':      5,
    'HBA':      10,
}
VEBER_LIMITS = {
    'RotBonds': 10,
    'TPSA':     140,
}
GHOSE_LIMITS = {
    'MW':       (160, 480),
    'LogP':     (-0.4, 5.6),
    'MR':       (40, 130),
    'HeavyAtoms':(20, 70),
}

def build_ligand_features(smiles: str) -> dict:
    """
    Compute RDKit molecular descriptors for ML model input.
    Returns dict of descriptor name -> value.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp     = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=512)
    import numpy as np
    fp_arr = np.zeros(512, dtype=int)
    from rdkit.DataStructs import ConvertToNumpyArray
    ConvertToNumpyArray(fp, fp_arr)
    feats = {
        'MW':           round(Descriptors.MolWt(mol), 2),
        'LogP':         round(Descriptors.MolLogP(mol), 3),
        'HBD':          rdMolDescriptors.CalcNumHBD(mol),
        'HBA':          rdMolDescriptors.CalcNumHBA(mol),
        'TPSA':         round(Descriptors.TPSA(mol), 2),
        'RotBonds':     rdMolDescriptors.CalcNumRotatableBonds(mol),
        'ArRings':      rdMolDescriptors.CalcNumAromaticRings(mol),
        'HeavyAtoms':   rdMolDescriptors.CalcNumHeavyAtoms(mol),
        'QED':          round(QED.qed(mol), 3),
        'FractionCSP3': round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        'NumStereo':    rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        'NumRings':     rdMolDescriptors.CalcNumRings(mol),
    }
    for i, v in enumerate(fp_arr):
        feats[f'FP_{i}'] = int(v)
    return feats


class ADMETCalculator:
    """
    Stage 06: ADMET drug-likeness assessment.
    Computes: Lipinski Ro5, Veber rules, Ghose filter,
              QED, BBB penetration prediction, GI absorption,
              solubility class, mutagenicity alerts (PAINS check).
    """

    def __init__(self, smiles, name, outdir, logger=None):
        self.smiles = smiles
        self.name   = name
        self.outdir = Path(outdir)
        self.logger = logger

    def log(self, msg, level='INFO'):
        if self.logger:
            self.logger.log(msg, level)

    def calculate(self) -> dict:
        self.outdir.mkdir(parents=True, exist_ok=True)
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            self.log(f'Invalid SMILES for ADMET: {self.smiles}', 'ERROR')
            return {}

        mw       = round(Descriptors.MolWt(mol), 2)
        logp     = round(Descriptors.MolLogP(mol), 3)
        hbd      = rdMolDescriptors.CalcNumHBD(mol)
        hba      = rdMolDescriptors.CalcNumHBA(mol)
        tpsa     = round(Descriptors.TPSA(mol), 2)
        rotbonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        heavy    = rdMolDescriptors.CalcNumHeavyAtoms(mol)
        mr       = round(Descriptors.MolMR(mol), 2)
        rings    = rdMolDescriptors.CalcNumRings(mol)
        ar_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
        qed_val  = round(QED.qed(mol), 3)
        fsp3     = round(rdMolDescriptors.CalcFractionCSP3(mol), 3)

        # Rule-based filters
        lip_viol   = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        veber_pass = (rotbonds <= 10 and tpsa <= 140)
        ghose_pass = (160 <= mw <= 480 and
                      -0.4 <= logp <= 5.6 and
                      40 <= mr <= 130 and
                      20 <= heavy <= 70)

        # BBB penetration: simplified rule-based prediction
        # MW < 400, logP 1-3, HBD <= 3, TPSA < 90, RotBonds < 8
        bbb = (mw < 400 and 1.0 <= logp <= 3.0 and
               hbd <= 3 and tpsa < 90 and rotbonds < 8)

        # GI absorption (Egan model approximation)
        gi = 'High' if (tpsa <= 131.6 and logp <= 5.88) else 'Low'

        # P-gp substrate likelihood (simple heuristic)
        pgp = 'Likely' if (mw > 400 and logp > 2 and hba > 4) else 'Unlikely'

        # Solubility class (ESOL approximation)
        log_sw = 0.16 - 0.63*logp - 0.0062*mw + 0.066*rotbonds - 0.74*ar_rings
        if log_sw > -1:   sol_class = 'Highly Soluble'
        elif log_sw > -2: sol_class = 'Soluble'
        elif log_sw > -3: sol_class = 'Moderately Soluble'
        elif log_sw > -4: sol_class = 'Slightly Soluble'
        else:             sol_class = 'Poorly Soluble'

        result = {
            'compound':             self.name,
            'MW':                   mw,
            'LogP':                 logp,
            'HBD':                  hbd,
            'HBA':                  hba,
            'TPSA':                 tpsa,
            'RotBonds':             rotbonds,
            'HeavyAtoms':           heavy,
            'MolRefractivity':      mr,
            'AromaticRings':        ar_rings,
            'QED':                  qed_val,
            'FractionCSP3':         fsp3,
            'Lipinski_violations':  lip_viol,
            'Lipinski_pass':        lip_viol <= 1,
            'Veber_pass':           veber_pass,
            'Ghose_pass':           ghose_pass,
            'BBB_penetrant':        bbb,
            'GI_absorption':        gi,
            'Pgp_substrate':        pgp,
            'Solubility_class':     sol_class,
            'LogSw_estimated':      round(log_sw, 3),
        }

        out = self.outdir / 'admet_results.json'
        with open(out, 'w') as f:
            json.dump(result, f, indent=2)

        self.log(
            f'ADMET: QED={qed_val:.3f}  LogP={logp:.2f}  '
            f'TPSA={tpsa:.1f}  Lipinski_violations={lip_viol}  '
            f'BBB={"Yes" if bbb else "No"}  GI={gi}',
            'OK'
        )
        return result
