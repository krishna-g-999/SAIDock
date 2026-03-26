import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from saidock.utils.pubchem_client import name_to_smiles

class LigandPrep:
    """
    Stage 03: Ligand preparation.
    Accepts drug name (PubChem lookup), SMILES string, or SDF file.
    Generates 3D conformer (ETKDG + MMFF94), Gasteiger charges,
    and converts to PDBQT for AutoDock Vina.
    """

    def __init__(self, outdir, logger=None,
                 name=None, smiles=None, sdf_path=None):
        self.outdir   = Path(outdir)
        self.logger   = logger
        self._name    = name
        self._smiles  = smiles
        self._sdf     = sdf_path
        self.name     = None
        self.smiles   = None
        self.mw       = 0.0
        self._mol     = None
        self._sdf_out = None
        self._pdbqt   = None

    def log(self, msg, level='INFO'):
        if self.logger:
            self.logger.log(msg, level)

    def prepare(self) -> Path:
        self.outdir.mkdir(parents=True, exist_ok=True)
        self._resolve_input()
        self._build_3d()
        pdbqt = self._convert_to_pdbqt()
        self.log(
            f'Ligand ready: {self.name}  '
            f'MW={self.mw:.2f}  HBD={rdMolDescriptors.CalcNumHBD(self._mol)}  '
            f'HBA={rdMolDescriptors.CalcNumHBA(self._mol)}',
            'OK'
        )
        return pdbqt

    def _resolve_input(self):
        if self._smiles:
            self.smiles = self._smiles
            self.name   = self._name or 'Compound'
        elif self._name:
            self.log(f'Fetching SMILES from PubChem: {self._name}')
            result      = name_to_smiles(self._name)
            self.smiles = result['smiles']
            self.name   = result['iupac_name'] or self._name
            self.log(f'PubChem CID: {result["cid"]}  '
                     f'SMILES: {self.smiles[:60]}')
        elif self._sdf:
            supplier = Chem.SDMolSupplier(str(self._sdf), removeHs=False)
            mol = next(iter(supplier), None)
            if mol is None:
                raise ValueError(f'Cannot read SDF: {self._sdf}')
            self._mol   = Chem.RemoveHs(mol)
            self.smiles = Chem.MolToSmiles(self._mol)
            self.name   = (mol.GetProp('_Name')
                           if mol.HasProp('_Name') else Path(self._sdf).stem)
        else:
            raise ValueError('Provide name, smiles, or sdf_path')

    def _build_3d(self):
        if self._mol is None:
            mol = Chem.MolFromSmiles(self.smiles)
            if mol is None:
                raise ValueError(f'Invalid SMILES: {self.smiles}')
            mol = Chem.AddHs(mol)
            ps  = AllChem.ETKDGv3()
            ps.randomSeed = 42
            result = AllChem.EmbedMolecule(mol, ps)
            if result == -1:
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
            self._mol = Chem.RemoveHs(mol)
        self.mw = Descriptors.MolWt(self._mol)
        sdf_path = self.outdir / f'{self.name}_3D.sdf'
        writer = Chem.SDWriter(str(sdf_path))
        writer.write(self._mol)
        writer.close()
        self._sdf_out = sdf_path
        self.log(f'3D conformer generated (ETKDG + MMFF94): {sdf_path.name}')

    def _convert_to_pdbqt(self) -> Path:
        out = self.outdir / f'{self.name}_ligand.pdbqt'
        ret = subprocess.run([
            'obabel', str(self._sdf_out), '-O', str(out),
            '--partialcharge', 'gasteiger', '-h',
            '--errorlevel', '1',
        ], capture_output=True, text=True)
        if out.exists() and out.stat().st_size > 50:
            self._pdbqt = out
            self.log(f'Ligand PDBQT generated: {out.name}', 'OK')
            return out
        raise FileNotFoundError(
            f'Ligand PDBQT conversion failed. '
            f'Check obabel installation.\n{ret.stderr[:200]}'
        )
