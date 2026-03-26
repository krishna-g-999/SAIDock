import json
import re
import subprocess
import numpy as np
from pathlib import Path


class MultiPocketDocker:
    """
    Stage 04: Multi-pocket molecular docking.
    Primary:  AutoDock Vina 1.2+
    Fallback: PolarDock (pure Python)
    """

    def __init__(self, receptor_pdbqt, ligand_pdbqt, pockets,
                 outdir, exhaustiveness=32, cpu=4,
                 n_poses=9, logger=None):
        self.receptor_pdbqt = Path(receptor_pdbqt)
        self.ligand_pdbqt   = Path(ligand_pdbqt)
        self.pockets        = pockets
        self.outdir         = Path(outdir)
        self.exhaustiveness = exhaustiveness
        self.cpu            = cpu
        self.n_poses        = n_poses
        self.logger         = logger

    def log(self, msg, level='INFO'):
        if self.logger:
            self.logger.log(msg, level)

    def run(self) -> list:
        self.outdir.mkdir(parents=True, exist_ok=True)
        vina_available = self._check_vina()
        self.log(f'Docking engine: {"AutoDock Vina" if vina_available else "PolarDock (built-in)"}')

        results = []
        for pk in self.pockets:
            pk_id = pk.get('pocket_id', 0)
            self.log(
                f'Docking pocket {pk_id}  '
                f'center={[round(x,1) for x in pk["center"]]}  '
                f'box={[round(b,1) for b in pk.get("box_size",[22,22,22])]}'
            )
            res = self._run_vina(pk) if vina_available else self._run_polar_dock(pk)
            res['pocket_id']          = pk_id
            res['pocket_drugg_score'] = pk.get('composite_score', 0)
            res['center']             = pk['center']
            res['engine']             = 'AutoDock Vina' if vina_available else 'PolarDock'
            results.append(res)
            self.log(
                f'  Pocket {pk_id}: best={res["best_score"]:.3f} kcal/mol  '
                f'n_poses={res["n_poses"]}',
                'OK'
            )

        results.sort(key=lambda r: r['best_score'])
        with open(self.outdir / 'docking_results.json', 'w') as f:
            json.dump(results, f, indent=2)
        return results

    def _check_vina(self) -> bool:
        try:
            ret = subprocess.run(['vina', '--version'],
                                 capture_output=True, text=True, timeout=5)
            if ret.returncode == 0:
                self.log(f'Vina version: {ret.stdout.strip()[:40]}')
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
        return False

    def _run_vina(self, pocket: dict) -> dict:
        pk_id      = pocket.get('pocket_id', 0)
        cx, cy, cz = pocket['center']

        # Use sensible default box size if pocket didn't compute one
        raw_box = pocket.get('box_size', [22, 22, 22])
        if isinstance(raw_box, list) and len(raw_box) == 3:
            sx = max(15.0, min(float(raw_box[0]), 40.0))
            sy = max(15.0, min(float(raw_box[1]), 40.0))
            sz = max(15.0, min(float(raw_box[2]), 40.0))
        else:
            sx = sy = sz = 22.0

        out_pose = self.outdir / f'pocket{pk_id}_poses.pdbqt'
        out_log  = self.outdir / f'pocket{pk_id}_vina.log'

        cmd = [
            'vina',
            '--receptor',       str(self.receptor_pdbqt),
            '--ligand',         str(self.ligand_pdbqt),
            '--out',            str(out_pose),
            '--center_x',       f'{cx:.3f}',
            '--center_y',       f'{cy:.3f}',
            '--center_z',       f'{cz:.3f}',
            '--size_x',         f'{sx:.1f}',
            '--size_y',         f'{sy:.1f}',
            '--size_z',         f'{sz:.1f}',
            '--exhaustiveness', str(self.exhaustiveness),
            '--num_modes',      str(self.n_poses),
            '--cpu',            str(self.cpu),
        ]

        self.log(f'  Vina cmd: {" ".join(cmd[6:14])}')  # log box params
        ret = subprocess.run(cmd, capture_output=True, text=True, timeout=max(600, int(getattr(self,"_vina_timeout",600))))
        full_log = ret.stdout + '\n' + ret.stderr
        out_log.write_text(full_log)

        if ret.returncode != 0:
            self.log(
                f'  Vina returned code {ret.returncode}. '
                f'Stderr: {ret.stderr[:200]}', 'WARN'
            )

        # Parse scores — try all methods in order
        scores = self._parse_vina_log_v12(full_log)
        if not scores:
            scores = self._parse_vina_log_v11(full_log)
        if not scores and out_pose.exists():
            scores = self._parse_pdbqt_remarks(out_pose)
        if not scores:
            self.log(
                f'  Score parsing failed. Log preview:\n'
                f'  {full_log[:400]}', 'WARN'
            )

        return self._summarise(scores, pk_id, str(out_pose))

    def _parse_vina_log_v12(self, log_text: str) -> list:
        """
        AutoDock Vina 1.2.x output format:
          mode |   affinity | dist from best mode
               | (kcal/mol) | rmsd l.b.| rmsd u.b.
          -----+------------+----------+----------
             1         -8.2      0.000      0.000
        """
        scores = []
        in_table = False
        for line in log_text.splitlines():
            if '-----+------------' in line or 'mode |   affinity' in line:
                in_table = True
                continue
            if in_table:
                m = re.match(r'\s+(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)', line)
                if m:
                    try:
                        scores.append(float(m.group(2)))
                    except ValueError:
                        pass
                elif line.strip() == '' and scores:
                    break
        return scores

    def _parse_vina_log_v11(self, log_text: str) -> list:
        """
        AutoDock Vina 1.1.x / alternative format:
          1    -8.2   0.000   0.000
        """
        scores = []
        for line in log_text.splitlines():
            m = re.match(r'\s*(\d+)\s+([-][\d.]+)\s+[\d.]+\s+[\d.]+', line)
            if m:
                try:
                    v = float(m.group(2))
                    if -25.0 < v < 0:
                        scores.append(v)
                except ValueError:
                    pass
        return scores

    def _parse_pdbqt_remarks(self, pdbqt_path: Path) -> list:
        """
        Extract scores from PDBQT MODEL REMARK lines:
          REMARK VINA RESULT:   -8.2      0.000      0.000
        """
        scores = []
        with open(pdbqt_path) as f:
            for line in f:
                if 'VINA RESULT' in line:
                    parts = line.split()
                    for p in parts:
                        try:
                            v = float(p)
                            if -25.0 < v < 0:
                                scores.append(v)
                                break
                        except ValueError:
                            continue
        return scores

    def _run_polar_dock(self, pocket: dict) -> dict:
        from saidock.core.polar_dock import PolarDock

        receptor_pdb = str(self.receptor_pdbqt).replace('_receptor.pdbqt', '_fixed.pdb')
        if not Path(receptor_pdb).exists():
            receptor_pdb = str(self.receptor_pdbqt).replace('.pdbqt', '_fixed.pdb')
        if not Path(receptor_pdb).exists():
            receptor_pdb = str(self.receptor_pdbqt).replace('.pdbqt', '.pdb')

        smiles_file = self.outdir.parent / 'ligand' / 'ligand_smiles.txt'
        if not smiles_file.exists():
            raise FileNotFoundError('Ligand SMILES file not found for PolarDock')
        smiles = smiles_file.read_text().strip()

        box = pocket.get('box_size', [22, 22, 22])
        docker = PolarDock(
            receptor_pdb   = receptor_pdb,
            pocket_center  = pocket['center'],
            box_size       = box,
            exhaustiveness = min(self.exhaustiveness, 16),
            n_poses        = self.n_poses,
        )
        poses  = docker.dock(smiles)
        scores = [p['score'] for p in poses] if poses else []
        pk_id  = pocket.get('pocket_id', 0)
        return self._summarise(scores, pk_id, poses=poses)

    def _summarise(self, scores, pk_id, pose_file=None, poses=None) -> dict:
        # Filter positive scores — they indicate severe steric clashes
        valid = [s for s in (scores or []) if s < 0]
        if not valid:
            if scores and max(scores) > 0:
                self.log(
                    f"  Pocket {pk_id}: clash pose(s) discarded "
                    f"(max={max(scores):.1f}). No valid poses.", "WARN"
                )
            valid = [-4.0]   # sentinel — downstream DTSS treats this as non-binder
        s    = sorted(valid)
        top3 = s[:3]
        return {
            "best_score": round(s[0], 3),
            "mean_top3":  round(float(np.mean(top3)), 3),
            "all_scores": [round(x, 3) for x in s],
            "n_poses":    len(valid),
            "pose_file":  str(pose_file) if pose_file else None,
            "pocket_id":  pk_id,
            "has_clash":  any(sc >= 0 for sc in (scores or [])),
        }

# ---- patch: filter positive scores at bottom of _summarise ----
def _summarise_patched(self, scores, pk_id, pose_file=None, poses=None):
    # Positive scores = severe steric clash — discard completely
    valid = [s for s in scores if s < 0]
    if not valid:
        self.log(
            f'  Pocket {pk_id}: all poses have positive score '
            f'(min={min(scores):.2f}) — steric clash. Skipping pocket.', 'WARN'
        )
        valid = [-4.0]   # fallback sentinel
    s    = sorted(valid)
    top3 = s[:3]
    return {
        'best_score': round(s[0], 3),
        'mean_top3':  round(float(sum(top3)/len(top3)), 3),
        'all_scores': [round(x, 3) for x in s],
        'n_poses':    len(valid),
        'pose_file':  str(pose_file) if pose_file else None,
        'pocket_id':  pk_id,
        'has_clash':  len(valid) != len(scores),
    }

MultiPocketDocker._summarise = _summarise_patched

# ---- patch: filter positive scores at bottom of _summarise ----
def _summarise_patched(self, scores, pk_id, pose_file=None, poses=None):
    # Positive scores = severe steric clash — discard completely
    valid = [s for s in scores if s < 0]
    if not valid:
        self.log(
            f'  Pocket {pk_id}: all poses have positive score '
            f'(min={min(scores):.2f}) — steric clash. Skipping pocket.', 'WARN'
        )
        valid = [-4.0]   # fallback sentinel
    s    = sorted(valid)
    top3 = s[:3]
    return {
        'best_score': round(s[0], 3),
        'mean_top3':  round(float(sum(top3)/len(top3)), 3),
        'all_scores': [round(x, 3) for x in s],
        'n_poses':    len(valid),
        'pose_file':  str(pose_file) if pose_file else None,
        'pocket_id':  pk_id,
        'has_clash':  len(valid) != len(scores),
    }

MultiPocketDocker._summarise = _summarise_patched
