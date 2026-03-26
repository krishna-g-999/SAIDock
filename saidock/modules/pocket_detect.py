import json
import subprocess
import numpy as np
from pathlib import Path
from saidock.core.polar_pocket import PolarPocket, read_protein_atoms

# Volume limits for a druggable pocket
# < 100 A3  : too small for any drug-like molecule
# > 10000 A3: too large — likely the entire protein surface, not a true pocket
POCKET_VOL_MIN  = 100.0
POCKET_VOL_MAX  = 10000.0

# Docking box size limits per dimension (Angstroms)
BOX_MIN = 18.0
BOX_MAX = 35.0


class PocketDetector:
    def __init__(self, receptor_pdb, outdir, n_pockets=5, logger=None):
        self.receptor_pdb = Path(receptor_pdb)
        self.outdir       = Path(outdir)
        self.n_pockets    = n_pockets
        self.logger       = logger

    def log(self, msg, level='INFO'):
        if self.logger:
            self.logger.log(msg, level)

    def detect(self) -> list:
        self.outdir.mkdir(parents=True, exist_ok=True)
        pockets = self._try_fpocket()
        if not pockets:
            self.log('fpocket not found or failed. Running PolarPocket (built-in engine).', 'WARN')
            pockets = self._run_polar_pocket()
        if not pockets:
            self.log('No pockets found. Using blind docking box.', 'WARN')
            pockets = [self._blind_box()]

        # Filter out unreasonable pockets (entire protein surface artefacts)
        valid = [p for p in pockets
                 if POCKET_VOL_MIN <= p.get('volume', 500) <= POCKET_VOL_MAX]
        if not valid:
            self.log(
                f'All {len(pockets)} pockets failed volume filter '
                f'({POCKET_VOL_MIN}-{POCKET_VOL_MAX} A3). '
                f'Using top 5 unfiltered by druggability.', 'WARN'
            )
            valid = pockets[:5]

        ranked = self._rank(valid)

        # Enforce sensible box sizes on every pocket
        for pk in ranked:
            pk['box_size'] = self._clamp_box(pk.get('box_size', [22, 22, 22]))

        out = self.outdir / 'pockets.json'
        with open(out, 'w') as f:
            json.dump(ranked, f, indent=2)

        for i, pk in enumerate(ranked[:3]):
            self.log(
                f'Pocket {i+1}: center={[round(x,1) for x in pk["center"]]}  '
                f'vol={pk.get("volume",0):.0f} A3  '
                f'box={[round(b,1) for b in pk["box_size"]]}  '
                f'drugg={pk.get("druggability_score",0):.3f}',
                'INFO'
            )
        return ranked[:self.n_pockets]

    @staticmethod
    def _clamp_box(box_size: list) -> list:
        """
        Compute a sensible per-pocket box size.
        Rules:
          - Minimum 20A (fits any drug-like molecule)
          - Maximum 35A (larger than any realistic binding site)
          - Use pocket dimension + 6A padding if PolarPocket computed it
          - If box is suspiciously small (< 20A) or large (> 35A), apply limits
        """
        BOX_HARD_MIN = 20.0
        BOX_HARD_MAX = 35.0
        result = []
        for b in box_size:
            b = float(b)
            # If PolarPocket returned a box from actual sphere span
            padded = b + 6.0
            clamped = max(BOX_HARD_MIN, min(padded, BOX_HARD_MAX))
            result.append(round(clamped, 1))
        return result

    def _try_fpocket(self) -> list:
        """
        Run fpocket on the fixed PDB receptor.

        fpocket always writes its output NEXT TO the input file, not in cwd.
        Strategy: copy the fixed PDB into self.outdir so that fpocket output
        lands in self.outdir/{stem}_out/ — exactly where the parser expects it.
        """
        import shutil

        # Copy fixed PDB into pocket outdir so fpocket output lands here
        local_pdb = self.outdir / self.receptor_pdb.name
        if not local_pdb.exists():
            shutil.copy2(str(self.receptor_pdb), str(local_pdb))

        out_dir   = self.outdir / f'{local_pdb.stem}_out'
        info_file = out_dir    / f'{local_pdb.stem}_info.txt'

        # Only run fpocket if output doesn't already exist
        if not info_file.exists():
            ret = subprocess.run(
                ['fpocket', '-f', str(local_pdb)],
                capture_output=True, text=True,
                cwd=str(self.outdir)
            )
            if ret.returncode != 0:
                self.log(f'fpocket error: {ret.stderr[:200]}', 'WARN')
                return []

        if not info_file.exists():
            return []

        pockets = self._parse_fpocket_output(info_file, out_dir)
        self.log(f'fpocket detected {len(pockets)} pockets', 'OK')
        return pockets

    def _parse_fpocket_output(self, info_file, out_dir) -> list:
        import re
        pockets, current = [], {}
        prop_map = {
            'druggability score':             'druggability_score',
            'volume :':                       'volume',        # exact match avoids 'volume score:'
            'mean local hydrophobic density': 'hydrophobicity',
            'number of alpha spheres':        'n_spheres',
            'total sasa':                     'total_sasa',
            'apolar sasa':                    'apolar_sasa',
        }
        with open(info_file) as f:
            for line in f:
                m = re.match(r'Pocket (\d+)\s*:', line)
                if m:
                    if current:
                        pockets.append(current)
                    current = {'pocket_id': int(m.group(1)), 'engine': 'fpocket'}
                for key, attr in prop_map.items():
                    if key in line.lower():
                        nums = re.findall(r'[-+]?\d*\.?\d+', line.split(':', 1)[-1])
                        if nums:
                            try:
                                current[attr] = float(nums[0])
                            except ValueError:
                                pass
        if current:
            pockets.append(current)
        for pk in pockets:
            pk_file = out_dir / 'pockets' / f'pocket{pk["pocket_id"]}_atm.pdb'
            if pk_file.exists():
                coords = self._coords_from_pdb(pk_file)
                if len(coords):
                    pk['center']   = coords.mean(axis=0).tolist()
                    rng = coords.max(axis=0) - coords.min(axis=0)
                    pk['box_size'] = (rng + 4).tolist()
        return pockets

    def _run_polar_pocket(self) -> list:
        coords, elements = read_protein_atoms(str(self.receptor_pdb))
        if len(coords) < 10:
            return []
        engine  = PolarPocket(coords, elements)
        pockets = engine.detect()
        for pk in pockets:
            pk['engine'] = 'PolarPocket'
        self.log(f'PolarPocket detected {len(pockets)} pockets')
        return pockets

    def _rank(self, pockets: list) -> list:
        for pk in pockets:
            ds    = pk.get('druggability_score', 0) or 0
            vol   = pk.get('volume', 300) or 300
            hydro = pk.get('hydrophobicity', 0.5) or 0.5
            n_sph = pk.get('n_spheres', 10) or 10
            vol_n   = min(vol / 1000.0, 1.0)
            hydro_n = min(hydro / 5.0, 1.0)
            sph_n   = min(n_sph / 30.0, 1.0)
            pk['composite_score'] = round(
                0.50*ds + 0.20*vol_n + 0.20*hydro_n + 0.10*sph_n, 4
            )
        return sorted(pockets,
                      key=lambda x: x.get('composite_score', 0),
                      reverse=True)

    def _blind_box(self) -> dict:
        coords = self._coords_from_pdb(self.receptor_pdb)
        if not len(coords):
            return {'pocket_id': 0, 'center': [0, 0, 0],
                    'box_size': [30, 30, 30], 'composite_score': 0,
                    'druggability_score': 0, 'volume': 500, 'engine': 'blind'}
        center   = coords.mean(axis=0).tolist()
        return {'pocket_id': 0, 'center': center,
                'box_size': [30, 30, 30],
                'composite_score': 0.2, 'druggability_score': 0.2,
                'volume': 500, 'engine': 'blind'}

    @staticmethod
    def _coords_from_pdb(pdb_path) -> np.ndarray:
        coords = []
        with open(pdb_path) as f:
            for line in f:
                if line[:4] in ('ATOM', 'HETA'):
                    try:
                        coords.append([float(line[30:38]),
                                       float(line[38:46]),
                                       float(line[46:54])])
                    except ValueError:
                        pass
        return np.array(coords) if coords else np.array([]).reshape(0, 3)
