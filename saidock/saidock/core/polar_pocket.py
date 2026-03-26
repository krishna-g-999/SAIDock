#!/usr/bin/env python3
"""
PolarPocket: Pure Python pocket detection engine.
Algorithm: Delaunay triangulation -> alpha sphere extraction
           -> single-linkage clustering -> druggability scoring
No external binaries required.
"""
import numpy as np
from scipy.spatial import Delaunay, cKDTree
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from typing import List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

ALPHA_MIN_RADIUS  = 3.0
ALPHA_MAX_RADIUS  = 7.5
CLUSTER_CUTOFF    = 4.5
MIN_CLUSTER_SIZE  = 5
VDW_RADII = {
    'C':1.70,'N':1.55,'O':1.52,'S':1.80,'P':1.80,'H':1.20,
    'F':1.47,'CL':1.75,'BR':1.85,'I':1.98,'ZN':1.39,
    'MG':1.73,'CA':1.97,'FE':1.60,'DEFAULT':1.70,
}
HYDROPHOBIC_ATOMS = {'C','S','F','CL','BR','I'}
POLAR_ATOMS       = {'N','O','S'}


class PolarPocket:
    def __init__(self, coords, elements,
                 min_radius=ALPHA_MIN_RADIUS,
                 max_radius=ALPHA_MAX_RADIUS,
                 cluster_cutoff=CLUSTER_CUTOFF,
                 min_cluster_size=MIN_CLUSTER_SIZE):
        self.coords           = np.array(coords, dtype=np.float64)
        self.elements         = [e.upper().strip() for e in elements]
        self.min_radius       = min_radius
        self.max_radius       = max_radius
        self.cluster_cutoff   = cluster_cutoff
        self.min_cluster_size = min_cluster_size

    def detect(self) -> List[Dict]:
        alphas = self._compute_alpha_spheres()
        if len(alphas) == 0:
            return []
        labels  = self._cluster(alphas)
        pockets = self._build_descriptors(labels, alphas)
        pockets = self._score(pockets)
        return sorted(pockets, key=lambda p: p['composite_score'], reverse=True)

    def _compute_alpha_spheres(self) -> np.ndarray:
        if len(self.coords) < 4:
            return np.array([]).reshape(0, 4)
        try:
            tri = Delaunay(self.coords)
        except Exception:
            return np.array([]).reshape(0, 4)
        tree    = cKDTree(self.coords)
        results = []
        for simplex in tri.simplices:
            pts = self.coords[simplex]
            center, radius = self._circumsphere(pts)
            if center is None:
                continue
            if self.min_radius <= radius <= self.max_radius:
                nn_dist, _ = tree.query(center, k=1)
                if nn_dist >= radius - 0.01:
                    results.append([center[0], center[1], center[2], radius])
        return np.array(results) if results else np.array([]).reshape(0, 4)

    @staticmethod
    def _circumsphere(pts):
        A = 2 * (pts[1:] - pts[0])
        b = np.sum(pts[1:]**2 - pts[0]**2, axis=1)
        try:
            center = np.linalg.solve(A, b)
            radius = np.linalg.norm(center - pts[0])
            return center, radius
        except np.linalg.LinAlgError:
            return None, None

    def _cluster(self, alphas) -> np.ndarray:
        if len(alphas) < 2:
            return np.zeros(len(alphas), dtype=int)
        dists  = pdist(alphas[:, :3])
        Z      = linkage(dists, method='single')
        return fcluster(Z, t=self.cluster_cutoff, criterion='distance')

    def _build_descriptors(self, labels, alphas) -> List[Dict]:
        pockets = []
        for cid in np.unique(labels):
            mask    = labels == cid
            spheres = alphas[mask]
            if len(spheres) < self.min_cluster_size:
                continue
            centers  = spheres[:, :3]
            radii    = spheres[:, 3]
            centroid = centers.mean(axis=0)
            volume   = float(np.sum((4/3) * np.pi * radii**3))
            lining   = self._lining_atoms(centers)
            if len(lining):
                elems     = [self.elements[i] for i in lining]
                n_hydro   = sum(1 for e in elems if e in HYDROPHOBIC_ATOMS)
                hydrophob = n_hydro / len(elems)
                polar_fr  = sum(1 for e in elems if e in POLAR_ATOMS) / len(elems)
            else:
                hydrophob = 0.5
                polar_fr  = 0.3
            rng = centers.max(axis=0) - centers.min(axis=0)
            pockets.append({
                'pocket_id':      int(cid),
                'n_spheres':      int(len(spheres)),
                'center':         centroid.tolist(),
                'volume':         round(volume, 2),
                'radius_mean':    round(float(radii.mean()), 3),
                'hydrophobicity': round(hydrophob, 3),
                'polar_fraction': round(polar_fr, 3),
                'n_atoms_lining': len(lining),
                'box_size':       [round(float(r + 4), 1) for r in rng],
            })
        return pockets

    def _lining_atoms(self, sphere_centers, cutoff=5.0):
        tree   = cKDTree(self.coords)
        lining = set()
        for sc in sphere_centers:
            lining.update(tree.query_ball_point(sc, r=cutoff))
        return sorted(lining)

    def _score(self, pockets) -> List[Dict]:
        for pk in pockets:
            vol   = pk['volume']
            n_sph = pk['n_spheres']
            hydro = pk['hydrophobicity']
            vol_s   = 1 / (1 + np.exp(-(vol - 300) / 200))
            sph_s   = min(n_sph / 30.0, 1.0)
            bal     = 1.0 - abs(hydro - 0.5) * 2.0
            hydro_s = min(hydro / 0.6, 1.0)
            ds = (0.35 * vol_s + 0.30 * hydro_s
                + 0.20 * sph_s + 0.15 * bal)
            pk['druggability_score'] = round(ds, 4)
            pk['composite_score']   = round(ds, 4)
        return pockets


def read_protein_atoms(pdb_path: str, chain: str = None):
    coords, elements = [], []
    with open(pdb_path) as f:
        for line in f:
            if line[:4] != 'ATOM':
                continue
            if chain and line[21] != chain.upper():
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = line[76:78].strip()
                if not elem:
                    elem = line[12:16].strip().lstrip('0123456789')[:2]
                coords.append([x, y, z])
                elements.append(elem.upper())
            except (ValueError, IndexError):
                continue
    return np.array(coords), elements
