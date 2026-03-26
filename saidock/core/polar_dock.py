#!/usr/bin/env python3
"""
PolarDock: Pure Python docking engine.
Scoring function: Vina empirical terms (Trott & Olson, 2010)
Search: Iterated Local Search + BFGS
No external binaries required.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from typing import List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

VINA_WEIGHTS = {
    'gauss1':     -0.0356,
    'gauss2':     -0.00516,
    'repulsion':   0.840,
    'hydrophobic':-0.0351,
    'hbond':      -0.587,
    'torsion':     0.0585,
}
VDW_RADII = {
    'C':1.90,'N':1.80,'O':1.70,'S':2.00,'P':2.10,'H':1.00,
    'F':1.54,'CL':1.77,'BR':1.90,'I':2.05,'ZN':1.44,
    'MG':1.73,'CA':1.80,'FE':1.60,'DEFAULT':1.80,
}

def get_vdw(symbol: str) -> float:
    return VDW_RADII.get(symbol.upper(), VDW_RADII['DEFAULT'])

def f_gauss(d, offset, width):
    return np.exp(-((d - offset) / width) ** 2)

def f_repulsion(d):
    return np.where(d < 0, d**2, 0.0)

def f_hydrophobic(d):
    return np.where(d <= 0.5, 1.0,
           np.where(d >= 1.5, 0.0, 1.0 - (d - 0.5)))

def f_hbond(d):
    return np.where(d <= -0.7, 1.0,
           np.where(d >= 0.0,  0.0, d / (-0.7)))


class VinaScorer:
    def __init__(self, rec_coords, rec_elems,
                 rec_donor, rec_acceptor, rec_hydrophobic,
                 cutoff=8.0):
        self.rc  = rec_coords
        self.re  = rec_elems
        self.rd  = rec_donor
        self.ra  = rec_acceptor
        self.rh  = rec_hydrophobic
        self.rv  = np.array([get_vdw(e) for e in rec_elems])
        self.cut = cutoff

    def score(self, lig_c, lig_e, lig_d, lig_a, lig_h, n_rot=0):
        lig_v = np.array([get_vdw(e) for e in lig_e])
        g1 = g2 = rep = hydr = hb = 0.0
        for i, (lc, lv, lhp, ld, la) in enumerate(
                zip(lig_c, lig_v, lig_h, lig_d, lig_a)):
            diffs   = self.rc - lc
            dist_sq = (diffs**2).sum(axis=1)
            within  = dist_sq < self.cut**2
            if not within.any():
                continue
            dist = np.sqrt(dist_sq[within])
            d    = dist - (self.rv[within] + lv)
            g1  += f_gauss(d, 0.0, 0.5).sum()
            g2  += f_gauss(d, 3.0, 2.0).sum()
            rep += f_repulsion(d).sum()
            if lhp:
                mask = self.rh[within]
                if mask.any():
                    hydr += f_hydrophobic(d[mask]).sum()
            if ld:
                mask = self.ra[within]
                if mask.any():
                    hb += f_hbond(d[mask]).sum()
            if la:
                mask = self.rd[within]
                if mask.any():
                    hb += f_hbond(d[mask]).sum()
        w = VINA_WEIGHTS
        return round(
            w['gauss1']*g1 + w['gauss2']*g2 + w['repulsion']*rep
            + w['hydrophobic']*hydr + w['hbond']*hb
            + w['torsion']*n_rot, 4
        )


class PolarDock:
    def __init__(self, receptor_pdb, pocket_center, box_size,
                 exhaustiveness=8, n_poses=9, seed=42):
        self.receptor_pdb  = receptor_pdb
        self.center        = np.array(pocket_center, dtype=np.float64)
        self.box_half      = np.array(box_size, dtype=np.float64) / 2.0
        self.exhaustiveness= exhaustiveness
        self.n_poses       = n_poses
        self.rng           = np.random.default_rng(seed)
        self._scorer       = None
        self._load_receptor()

    def _load_receptor(self):
        from saidock.core.polar_pocket import read_protein_atoms
        coords, elems = read_protein_atoms(self.receptor_pdb)
        donor  = np.array([e in ('N','O') for e in elems])
        accept = np.array([e in ('N','O','F') for e in elems])
        hydro  = np.array([e in ('C','S') for e in elems])
        self._scorer = VinaScorer(coords, elems, donor, accept, hydro)

    def dock(self, smiles: str) -> List[Dict]:
        mol = self._prepare_ligand(smiles)
        if mol is None:
            return []
        lig_c0, le, ld, la, lh, n_rot = self._lig_features(mol)
        poses = []
        for run in range(self.exhaustiveness):
            x0  = self._random_start(lig_c0)
            pv  = np.array([*self.center, 0.0, 0.0, 0.0])
            bv, bs = self._ils(pv, lig_c0, le, ld, la, lh, n_rot)
            cf  = self._apply_pose(lig_c0, bv)
            poses.append({'score': bs, 'coords': cf})
        return self._cluster_poses(poses)[:self.n_poses]

    def _prepare_ligand(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        ps  = AllChem.ETKDGv3()
        ps.randomSeed = 42
        AllChem.EmbedMolecule(mol, ps)
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.RemoveHs(mol)

    def _lig_features(self, mol):
        conf  = mol.GetConformer()
        c     = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
        e     = [mol.GetAtomWithIdx(i).GetSymbol() for i in range(mol.GetNumAtoms())]
        d     = np.array([mol.GetAtomWithIdx(i).GetAtomicNum() in (7,8) and
                          mol.GetAtomWithIdx(i).GetTotalNumHs()>0
                          for i in range(mol.GetNumAtoms())])
        a     = np.array([mol.GetAtomWithIdx(i).GetAtomicNum() in (7,8,9)
                          for i in range(mol.GetNumAtoms())])
        h     = np.array([mol.GetAtomWithIdx(i).GetAtomicNum() == 6
                          for i in range(mol.GetNumAtoms())])
        n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        return c, e, d, a, h, n_rot

    def _random_start(self, c0):
        shift = (self.rng.random(3)*2 - 1) * self.box_half * 0.7
        return c0 - c0.mean(axis=0) + self.center + shift

    def _apply_pose(self, c0, pv):
        tx,ty,tz,rx,ry,rz = pv
        rot     = Rotation.from_rotvec([rx, ry, rz])
        cent    = c0.mean(axis=0)
        return rot.apply(c0 - cent) + cent + np.array([tx,ty,tz]) - cent

    def _ils(self, pv0, c0, le, ld, la, lh, n_rot, n_steps=15):
        def obj(pv):
            return self._scorer.score(self._apply_pose(c0, pv),
                                      le, ld, la, lh, n_rot)
        cur_pv = pv0.copy()
        cur_sc = obj(cur_pv)
        bv, bs = cur_pv.copy(), cur_sc
        for step in range(n_steps):
            res    = minimize(obj, cur_pv, method='BFGS',
                              options={'maxiter':30,'gtol':1.0})
            new_pv = res.x
            new_sc = res.fun
            kT     = 1.0 * (0.1/1.0) ** (step/n_steps)
            delta  = new_sc - cur_sc
            if delta < 0 or self.rng.random() < np.exp(-delta/kT):
                cur_pv, cur_sc = new_pv, new_sc
            if cur_sc < bs:
                bv, bs = cur_pv.copy(), cur_sc
            perturb  = np.concatenate([self.rng.normal(0,1.0,3),
                                       self.rng.normal(0,0.3,3)])
            cur_pv   = cur_pv + perturb
            cur_pv[:3] = np.clip(cur_pv[:3],
                                  self.center - self.box_half,
                                  self.center + self.box_half)
        return bv, bs

    def _cluster_poses(self, poses, rmsd_cut=2.0):
        sorted_p = sorted(poses, key=lambda p: p['score'])
        clusters = []
        for pose in sorted_p:
            new = True
            for cl in clusters:
                n = min(len(pose['coords']), len(cl['coords']))
                rmsd = float(np.sqrt(
                    ((pose['coords'][:n]-cl['coords'][:n])**2
                     ).sum(axis=1).mean()))
                if rmsd < rmsd_cut:
                    new = False
                    break
            if new:
                clusters.append(pose)
        for i, cl in enumerate(clusters):
            cl['rank']  = i+1
            cl['score'] = round(cl['score'], 3)
        return clusters
