"""
SAIDock RNA-Targeted Docking Module
====================================
Docks small molecules to RNA-binding protein interfaces,
specifically the nucleic acid recognition cleft.

Scientific basis:
  - TDP-43 RRM1 binds RNA via Trp113/Arg171 cleft (PMID 31241884)
  - Compounds competing with RNA binding are validated ALS leads
  - rDock + AutoDock Vina both support RNA receptor preparation

Mode: Protein structure with RNA binding site as docking box
      (not RNA itself — protein is receptor, RNA site is cavity)
"""

from pathlib import Path
import json, subprocess, shutil, tempfile
from typing import Optional

# ── RNA-binding site definitions per target ───────────────────────────
# Key residues forming the RNA recognition cleft (from crystal structures)
# Used to define docking box centroid + radius
RNA_BINDING_SITES = {
    'TDP43_RRM1': {
        'pdb_id':       '5HDN',
        'center':       (18.5, 42.3, 25.1),   # RNP-1/RNP-2 centroid
        'box_size':     (20, 20, 20),           # Angstroms
        'key_residues': [109, 110, 112, 113, 171, 172, 174, 195],
        'note':         'Trp113/Arg171 cleft — competes with (UG)n RNA binding',
        'reference':    'PMID:31241884',
    },
    'TDP43_RRM2': {
        'pdb_id':       '2N3Y',
        'center':       (22.1, 38.7, 19.4),
        'box_size':     (20, 20, 20),
        'key_residues': [192, 194, 213, 259, 261],
        'note':         'RRM2 domain — second RNA recognition surface',
        'reference':    'PMID:28652374',
    },
    'Keap1_PPI': {
        'pdb_id':       '4XMB',
        'center':       (35.2, 55.8, 12.6),
        'box_size':     (22, 22, 22),
        'key_residues': [363, 364, 365, 407, 408, 410, 462, 473, 508, 509, 514],
        'note':         'Keap1 Kelch domain — Nrf2 ETGE/DLG motif binding groove',
        'reference':    'PMID:24818526',
    },
}

class RNADockingMode:
    """
    RNA-competition docking: dock small molecules to the nucleic acid
    binding interface of RNA-binding proteins.

    This is distinct from standard active-site docking — the box is
    centered on the RNA recognition cleft, not the protein core.
    Compounds scoring well here are RNA-protein interaction disruptors.
    """

    def __init__(self, vina_bin: str = 'vina', log_fn=None):
        self.vina_bin = vina_bin
        self.log      = log_fn or (lambda msg, lvl='INFO': print(f"[{lvl}] {msg}"))

    def get_rna_site(self, target_id: str) -> Optional[dict]:
        """Return RNA binding site definition for target."""
        # Try direct lookup
        site = RNA_BINDING_SITES.get(target_id)
        if site: return site
        # Fuzzy match
        tid = target_id.upper()
        for k, v in RNA_BINDING_SITES.items():
            if k.upper() in tid or tid in k.upper():
                return v
            if v.get('pdb_id','').upper() == tid:
                return v
        return None

    def dock_rna_competition(
        self,
        smiles:    str,
        target_id: str,
        receptor_pdbqt: str,
        ligand_pdbqt:   str,
        n_poses:   int = 9,
    ) -> dict:
        """
        Run AutoDock Vina with RNA-binding site box.
        Returns dict with rna_vina_score + rna_displacement_score.
        """
        site = self.get_rna_site(target_id)
        if site is None:
            self.log(f"No RNA site defined for {target_id}", 'WARN')
            return {'rna_mode': False, 'rna_vina_score': None}

        cx, cy, cz = site['center']
        sx, sy, sz = site['box_size']

        with tempfile.TemporaryDirectory() as tmpdir:
            out_pdbqt = Path(tmpdir) / 'rna_poses.pdbqt'
            vina_log  = Path(tmpdir) / 'vina_rna.log'

            cmd = [
                self.vina_bin,
                '--receptor', receptor_pdbqt,
                '--ligand',   ligand_pdbqt,
                '--center_x', str(cx),
                '--center_y', str(cy),
                '--center_z', str(cz),
                '--size_x',   str(sx),
                '--size_y',   str(sy),
                '--size_z',   str(sz),
                '--num_modes', str(n_poses),
                '--exhaustiveness', '16',
                '--out',  str(out_pdbqt),
                '--log',  str(vina_log),
            ]

            try:
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=300
                )
                scores = self._parse_vina_log(vina_log)
                best   = min(scores) if scores else None

                # RNA displacement score: how well does this molecule
                # compete with RNA? Compare to standard active site score.
                rna_score = self._compute_rna_displacement_score(
                    best, smiles, site
                )

                self.log(
                    f"RNA-mode docking [{target_id}]: "
                    f"best_score={best}  rna_disp={rna_score:.4f}", 'INFO'
                )
                return {
                    'rna_mode':               True,
                    'rna_binding_site':       site['note'],
                    'rna_key_residues':       site['key_residues'],
                    'rna_vina_score':         best,
                    'rna_displacement_score': rna_score,
                    'rna_all_scores':         scores[:3],
                    'reference':              site['reference'],
                }

            except subprocess.TimeoutExpired:
                self.log("RNA docking timed out", 'WARN')
                return {'rna_mode': False, 'rna_vina_score': None}
            except Exception as e:
                self.log(f"RNA docking error: {e}", 'WARN')
                return {'rna_mode': False, 'rna_vina_score': None}

    def _parse_vina_log(self, log_path: Path) -> list:
        scores = []
        if not log_path.exists(): return scores
        for line in log_path.read_text().splitlines():
            parts = line.strip().split()
            if parts and parts[0].isdigit():
                try: scores.append(float(parts[1]))
                except (ValueError, IndexError): pass
        return scores

    def _compute_rna_displacement_score(
        self,
        vina_score: Optional[float],
        smiles: str,
        site: dict,
    ) -> float:
        """
        RNA Displacement Score (RDS): probability that this compound
        would displace RNA from the binding cleft.

        Formula:
          RDS = 0.60 * DG_norm + 0.25 * shape_complementarity + 0.15 * charge_score

        DG_norm: normalised from Vina score (−12 kcal/mol = 1.0)
        Shape:   estimated from MW and planarity (aromatic rings)
        Charge:  positive charge favours RNA (phosphate) displacement
        """
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        import math

        # DG_norm
        if vina_score is None:
            dg_norm = 0.5
        else:
            dg_norm = min(1.0, max(0.0, abs(vina_score) / 12.0))

        # Shape complementarity proxy (planar aromatic molecules fit RRM cleft)
        try:
            mol        = Chem.MolFromSmiles(smiles)
            n_arom     = rdMolDescriptors.CalcNumAromaticRings(mol) if mol else 0
            mw         = Descriptors.MolWt(mol) if mol else 300
            shape_sc   = min(1.0, n_arom / 3.0) * min(1.0, mw / 400.0)
        except Exception:
            shape_sc = 0.5

        # Charge proxy (basic N atoms = positive charge = RNA competition)
        try:
            n_basic = sum(1 for a in mol.GetAtoms()
                         if a.GetAtomicNum()==7 and
                         a.GetTotalValence() < 4) if mol else 0
            charge_sc = min(1.0, n_basic / 3.0)
        except Exception:
            charge_sc = 0.3

        rds = 0.60 * dg_norm + 0.25 * shape_sc + 0.15 * charge_sc
        return round(min(1.0, max(0.0, rds)), 4)

    def score_rna_displacement_smiles_only(
        self, smiles: str, target_id: str
    ) -> dict:
        """
        Lightweight RNA displacement scoring without running Vina —
        uses only SMILES features + known RNA-binding pharmacophore rules.
        Use this when full docking is unavailable.
        """
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
        import numpy as np

        site = self.get_rna_site(target_id)
        if site is None:
            return {'rna_mode': False,
                    'rna_displacement_score': 0.5,
                    'rna_pharmacophore_match': False}

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'rna_mode': False,
                    'rna_displacement_score': 0.5,
                    'rna_pharmacophore_match': False}

        # RNA recognition cleft pharmacophore for TDP-43 RRM1:
        #   - Planar aromatic system (stacks with Trp113)
        #   - H-bond donor (contacts Arg171 backbone)
        #   - MW 200–500 (fits cleft without steric clash)
        n_arom  = rdMolDescriptors.CalcNumAromaticRings(mol)
        n_hbd   = rdMolDescriptors.CalcNumHBD(mol)
        n_hba   = rdMolDescriptors.CalcNumHBA(mol)
        mw      = Descriptors.MolWt(mol)
        tpsa    = rdMolDescriptors.CalcTPSA(mol)
        qed_val = Descriptors.qed(mol)

        # Pharmacophore match (simplified)
        pharm_match = (n_arom >= 1 and n_hbd >= 1 and 150 < mw < 550)

        # Weighted score
        arom_score  = min(1.0, n_arom / 3.0)
        hbd_score   = min(1.0, n_hbd  / 3.0)
        mw_score    = 1.0 - abs(mw - 320) / 320.0
        mw_score    = max(0.0, min(1.0, mw_score))
        tpsa_score  = max(0.0, 1.0 - tpsa / 140.0)

        rds_ligand = (0.35 * arom_score +
                      0.25 * hbd_score  +
                      0.25 * mw_score   +
                      0.15 * tpsa_score)

        return {
            'rna_mode':                True,
            'rna_binding_site':        site['note'],
            'rna_displacement_score':  round(rds_ligand, 4),
            'rna_pharmacophore_match': pharm_match,
            'rna_components': {
                'aromatic_score':   round(arom_score, 4),
                'hbd_score':        round(hbd_score,  4),
                'mw_score':         round(mw_score,   4),
                'tpsa_score':       round(tpsa_score, 4),
            },
            'reference': site['reference'],
        }
