#!/usr/bin/env python3
"""
SAIDock SurfaceAnalyser — Pan-pocket profiling with biological role classification.
Verified PolarPocket API: PolarPocket(pdb_path_str).detect()  → list of dicts
Pocket keys: pocket_id, center, volume, druggability_score, hydrophobicity,
             polar_fraction, n_spheres, box_size, composite_score
"""
import json
import math
import re
import traceback
import requests
from collections import Counter
from pathlib import Path

import numpy as np

_NGL_HEAD_SCRIPT = (
    "<script src='https://cdn.jsdelivr.net/npm/"
    "ngl@2.0.0-dev.37/dist/ngl.js'></script>"
)

# NGL body template — placeholders replaced at runtime
# Uses %s-style substitution, never f-string, to avoid brace conflicts
_NGL_BODY_TPL = (
    "<h2 style='color:#333;margin-top:20px'>3D Pocket Viewer</h2>"
    "<div id='sai-ngl' style='width:100%;height:460px;"
    "background:#0d1117;border:2px solid #30363d;"
    "border-radius:8px;margin:10px 0'></div>"
    "<div id='sai-ctrl' style='margin-bottom:10px'></div>"
    "<script>\n"
    "(function(){\n"
    "  var stage=new NGL.Stage('sai-ngl',{backgroundColor:'#0d1117'});\n"
    "  var pdat=__POCKETJSON__;\n"
    "  var recf='__RECPATH__';\n"
    "  if(recf){stage.loadFile(recf,{ext:'pdb'}).then(function(o){\n"
    "    o.addRepresentation('cartoon',{color:'chainid',opacity:0.85});\n"
    "    o.addRepresentation('surface',{opacity:0.1,color:'skyblue'});\n"
    "    stage.autoView();});}\n"
    "  var sh=[];\n"
    "  pdat.forEach(function(p,i){\n"
    "    var col=p.conf.indexOf('CONFIRM')>=0?[1,.2,.2]:\n"
    "            p.conf.indexOf('METAL')>=0?[.7,0,1]:\n"
    "            p.conf.indexOf('FUNCT')>=0?[1,.6,0]:[.2,.6,1];\n"
    "    var s=new NGL.Shape('p'+i);\n"
    "    s.addSphere([p.cx,p.cy,p.cz],col,3.2,'Pkt '+p.id);\n"
    "    var c=stage.addComponentFromObject(s);\n"
    "    c.addRepresentation('buffer');\n"
    "    sh.push({c:c,on:true});});\n"
    "  var ctrl=document.getElementById('sai-ctrl');\n"
    "  pdat.forEach(function(p,i){\n"
    "    var b=document.createElement('button');\n"
    "    b.style.cssText='background:#238636;border:1px solid #3fb950;"
    "color:#fff;padding:4px 10px;margin:3px;cursor:pointer;"
    "border-radius:4px;font-size:11px';\n"
    "    b.textContent='Pkt '+p.id+' ('+p.dtss+')';\n"
    "    b.onclick=(function(i,b){return function(){\n"
    "      sh[i].on=!sh[i].on;sh[i].c.setVisibility(sh[i].on);\n"
    "      b.style.background=sh[i].on?'#238636':'#21262d';\n"
    "    };})(i,b);ctrl.appendChild(b);});\n"
    "})();\n"
    "</script>"
)

# Amino acid property sets
POLAR_AA    = {'SER','THR','CYS','TYR','ASN','GLN','HIS'}
CHARGED_AA  = {'ARG','LYS','ASP','GLU'}
HYDRO_AA    = {'ALA','VAL','LEU','ILE','MET','PHE','TRP','PRO','GLY'}
AROMATIC_AA = {'PHE','TYR','TRP','HIS'}
METAL_AA    = {'HIS','CYS','ASP','GLU'}


class SurfaceAnalyser:
    _FPOCKET_KEY_MAP = {
        'volume'                        : 'volume',
        'pocket volume'                 : 'volume',
        'score'                         : 'score',
        'druggability score'            : 'druggability',
        'number of alpha spheres'       : 'n_alpha_spheres',
        'total sasa'                    : 'sasa_total',
        'polar sasa'                    : 'sasa_polar',
        'apolar sasa'                   : 'sasa_apolar',
        'hydrophobicity score'          : 'hydrophobicity',
        'volume score'                  : 'volume_score',
        'polarity score'                : 'polarity_score',
        'charge score'                  : 'charge_score',
        'flexibility'                   : 'flexibility',
        'alpha sphere density'          : 'alpha_density',
        'proportion of polar atoms'     : 'polar_proportion',
    }

    """Pan-pocket surface profiler for any PDB target."""

    MIN_VOL_PRIMARY     = 200.0
    MAX_VOL_ION         = 250.0
    MIN_VOL_PPI         = 700.0
    MIN_DRUGG_PRIMARY   = 0.55
    MIN_DRUGG_DRUGGABLE = 0.35
    ALLOSTERIC_MIN_DIST = 12.0
    METAL_CLUSTER_MIN   = 2
    VOL_NORMALISE_CAP   = 30000.0   # PolarPocket sums sphere volumes — cap display

    def __init__(self, receptor_pdb: str, target_id: str = '',
                 uniprot_id: str = None, outdir: str = '.',
                 min_vol: float = 100.0, max_pockets: int = 20,
                 logger=None):
        self.receptor_pdb = Path(receptor_pdb)
        self.target_id    = target_id
        self.uniprot_id   = uniprot_id
        self.outdir       = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.min_vol      = min_vol
        self.max_pockets  = max_pockets
        self.logger       = logger
        self._residues    = []

    # ── Logging ───────────────────────────────────────────────────────────────
    def log(self, msg: str, level: str = 'INFO'):
        if self.logger:
            self.logger.log(msg, level)
        else:
            print(f'  [{level}] {msg}')

    # ── Main entry ────────────────────────────────────────────────────────────
    def analyse(self) -> list:
        """Run complete surface analysis. Returns list of annotated pocket dicts."""
        try:
            return self._run_analysis()
        except Exception as e:
            self.log(f'Surface analysis error: {e}', 'WARN')
            self.log(traceback.format_exc(), 'WARN')
            return []

    def _run_analysis(self) -> list:
        # 1. Parse residues for composition analysis
        self._residues = self._parse_pdb_residues()
        self.log(f'Parsed {len(self._residues)} residues')

        # 2. Detect pockets
        raw = self._detect_pockets()
        self.log(f'Raw pockets detected: {len(raw)}')

        # 3. Filter + rank
        pockets = [p for p in raw if p['volume'] >= self.min_vol]
        pockets = sorted(pockets, key=lambda p: p['druggability'], reverse=True)
        pockets = pockets[:self.max_pockets]
        self.log(f'Pockets after filter (vol≥{self.min_vol} Å³): {len(pockets)}')

        if not pockets:
            return []

        # 4. Annotate each pocket
        primary_center = pockets[0]['center']
        for i, p in enumerate(pockets):
            p['rank']               = i + 1
            p['residue_profile']    = self._pocket_residue_profile(p)
            p['dist_from_primary']  = (0.0 if i == 0
                                        else self._dist(p['center'], primary_center))
            p['pocket_type']        = self._classify(p, i)
            p['interpretation']     = self._interpret(p)
            p['druggability_label'] = self._drugg_label(p['druggability'])
            # Propagate residue numbers for annotator
            p['_lining_resnums']    = p['residue_profile'].get('residue_numbers', [])

        # 5. Multi-evidence annotation (L1-L5)
        try:
            from saidock.modules.pocket_annotator import PocketAnnotator
            # Find raw PDB (before chain extraction) — one dir up from fixed
            raw_candidates = list(self.receptor_pdb.parent.glob('*_raw.pdb')) +                              list(self.receptor_pdb.parent.glob(
                                 f'{self.target_id}_raw.pdb'))
            raw_pdb = raw_candidates[0] if raw_candidates else self.receptor_pdb

            annotator = PocketAnnotator(
                pdb_id       = self.target_id,
                raw_pdb_path = str(raw_pdb),
                fixed_pdb_path = str(self._resolve_pdb_path() or self.receptor_pdb),
                chain        = 'A',
                logger       = self.logger,
                cache_dir    = str(self.outdir / '.cache'),
            )
            pockets = annotator.annotate_all(pockets)
        except Exception as e:
            self.log(f'PocketAnnotator failed (non-critical): {e}', 'WARN')

        # 6. UniProt cross-reference
        if self.uniprot_id or self.target_id:
            self._fetch_uniprot_annotations(pockets)

        # 6. Save atlas
        atlas_path = self.outdir / 'pocket_atlas.json'
        with open(atlas_path, 'w') as fh:
            json.dump(pockets, fh, indent=2, default=str)

        csv_path = self.outdir / 'pocket_atlas.csv'
        self._write_csv(pockets, csv_path)

        self.log(f'Pocket atlas saved: {atlas_path}', 'OK')
        # ── Write HTML report ────────────────────────────
        try:
            self.write_report(pockets)
        except Exception as _wr_e:
            import traceback as _tb
            err_log = self.outdir / 'report_error.log'
            err_log.write_text(_tb.format_exc())
            self.log(f'write_report failed: {_wr_e}  → {err_log}', 'ERROR')
        # ──────────────────────────────────────────────────
        return pockets

    # ── Pocket detection ──────────────────────────────────────────────────────
    def _detect_pockets(self) -> list:
        """Try fpocket first, fall back to PolarPocket."""
        pockets = self._run_fpocket()
        if pockets:
            self.log(f'fpocket: {len(pockets)} pockets', 'OK')
            return pockets
        self.log('fpocket not available — using PolarPocket', 'WARN')
        return self._run_polarpocket()

    def _run_fpocket(self) -> list:
        """Run fpocket and return pocket list with center coordinates.
        
        fpocket writes:
          {stem}_out/{stem}_info.txt          ← scores
          {stem}_out/pockets/pocket{N}_atm.pdb ← alpha spheres (used for center)
        """
        import shutil as _sh, subprocess, tempfile, numpy as _np
        if not _sh.which('fpocket'):
            return []
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_pdb = Path(tmpdir) / self.receptor_pdb.name
                _sh.copy(str(self.receptor_pdb), str(tmp_pdb))
                ret = subprocess.run(
                    ['fpocket', '-f', str(tmp_pdb)],
                    capture_output=True, text=True, timeout=120, cwd=tmpdir,
                )
                if ret.returncode not in (0, 1):
                    return []

                out_dir    = Path(tmpdir) / (tmp_pdb.stem + '_out')
                info_files = list(out_dir.rglob('*_info.txt'))
                pdb_dir    = out_dir / 'pockets'

                if not info_files:
                    return []

                pockets = self._parse_fpocket_info(info_files[0])

                # Compute pocket centres from alpha-sphere PDB files
                for p in pockets:
                    pid  = p.get('pocket_id', 0)
                    # fpocket names: pocket1_atm.pdb, pocket1_vert.pdb
                    atm  = pdb_dir / f'pocket{pid}_atm.pdb'
                    vert = pdb_dir / f'pocket{pid}_vert.pdb'
                    coords = []
                    for pf in [atm, vert]:
                        if pf.exists():
                            for line in pf.read_text().splitlines():
                                if line.startswith(('ATOM', 'HETATM')):
                                    try:
                                        coords.append([float(line[30:38]),
                                                       float(line[38:46]),
                                                       float(line[46:54])])
                                    except ValueError:
                                        pass
                            if coords:
                                break
                    if coords:
                        arr = _np.array(coords)
                        p['center'] = arr.mean(axis=0).tolist()
                    else:
                        # Fallback: no PDB — skip this pocket
                        p['center'] = None

                pockets = [p for p in pockets if p.get('center')]
                self.log(f'fpocket: {len(pockets)} pockets with centres', 'INFO')
                return pockets

        except Exception as e:
            self.log(f'fpocket error: {e}', 'WARN')
            return []

    def _parse_fpocket_info(self, info_file: Path) -> list:
        pockets, cur = [], {}
        with open(info_file) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('Pocket'):
                    if cur:
                        pockets.append(cur)
                    try:
                        cur = {'pocket_id': int(line.split()[1].rstrip(':'))}
                    except (ValueError, IndexError):
                        cur = {'pocket_id': len(pockets) + 1}
                elif ':' in line:
                    k, _, v = line.partition(':')
                    k_norm   = k.strip().lower()
                    internal = self._FPOCKET_KEY_MAP.get(
                        k_norm, k_norm.replace(' ', '_').replace('.', ''))
                    try:
                        cur[internal] = float(v.strip())
                    except (ValueError, TypeError):
                        cur[internal] = v.strip()
        if cur:
            pockets.append(cur)
        return pockets

    def _run_polarpocket(self) -> list:
        """
        Verified working call (from diagnostic):
          coords_arr = np.array([[x,y,z],...], dtype=float32)
          elements   = ['C','N','O',...]
          PolarPocket(coords_arr, elements).detect()  →  list of pocket dicts
        """
        try:
            import numpy as np
            from saidock.core.polar_pocket import PolarPocket

            pdb_path = self._resolve_pdb_path()
            if pdb_path is None:
                self.log(f'Cannot find PDB file near {self.receptor_pdb}', 'WARN')
                return []

            # Parse coords + elements from PDB (required by PolarPocket API)
            coords, elements = [], []
            with open(pdb_path) as fh:
                for line in fh:
                    if line[:4] not in ('ATOM', 'HETA'):
                        continue
                    try:
                        x  = float(line[30:38])
                        y  = float(line[38:46])
                        z  = float(line[46:54])
                        el = line[76:78].strip()
                        if not el:
                            el = line[12:16].strip().lstrip('0123456789')[:1]
                        coords.append([x, y, z])
                        elements.append(el or 'C')
                    except (ValueError, IndexError):
                        continue

            if not coords:
                self.log('No atoms parsed from PDB', 'WARN')
                return []

            self.log(f'PolarPocket using: {pdb_path.name} ({len(coords)} atoms)')
            coords_arr = np.array(coords, dtype=np.float32)

            # ← EXACT working call confirmed by diagnostic
            pp  = PolarPocket(coords_arr, elements)
            raw = pp.detect()

            if not raw:
                self.log('PolarPocket returned 0 pockets', 'WARN')
                return []

            self.log(f'PolarPocket: {len(raw)} pockets detected', 'OK')
            return self._normalise_pp_pockets(raw)

        except Exception as e:
            self.log(f'PolarPocket failed: {e}', 'WARN')
            self.log(traceback.format_exc().splitlines()[-1], 'WARN')
            return []

    def _resolve_pdb_path(self) -> Path:
        """Find a valid PDB file starting from self.receptor_pdb."""
        candidates = [self.receptor_pdb]
        # Search parent + grandparent for *_fixed.pdb and *_chainA.pdb
        for parent in [self.receptor_pdb.parent,
                       self.receptor_pdb.parent.parent]:
            candidates += list(parent.glob('*_fixed.pdb'))
            candidates += list(parent.glob('*_chainA.pdb'))
        for c in candidates:
            if Path(c).exists():
                return Path(c)
        return None

    def _resolve_chainA_pdb(self) -> Path:
        """
        Return chainA.pdb (pre-fix, author residue numbers preserved).
        This ensures residue numbering is consistent with B-factor proxy
        and ConSurfDB outputs, which use author numbering from the crystal.
        """
        for parent in [self.receptor_pdb.parent,
                       self.receptor_pdb.parent.parent]:
            for p in parent.glob('*_chainA.pdb'):
                return p
        return None

    @staticmethod
    def _normalise_pp_pockets(raw: list) -> list:
        """Map PolarPocket dict keys → standard SurfaceAnalyser schema."""
        VOL_CAP = 30000.0
        result  = []
        for i, r in enumerate(raw or []):
            if not isinstance(r, dict):
                continue
            center = r.get('center', [0.0, 0.0, 0.0])
            if not isinstance(center, (list, tuple)) or len(center) < 3:
                center = [0.0, 0.0, 0.0]

            vol = float(r.get('volume', 500.0))
            if vol > VOL_CAP:
                vol = min(vol / 100.0, VOL_CAP)

            drugg = float(
                r.get('druggability_score',
                r.get('composite_score',
                r.get('drug_score',
                r.get('druggability', 0.5))))
            )
            drugg = max(0.0, min(1.0, drugg))

            result.append({
                'pocket_id':      r.get('pocket_id', i + 1),
                'center':         [float(c) for c in center[:3]],
                'volume':         round(vol, 1),
                'druggability':   round(drugg, 4),
                'n_spheres':      r.get('n_spheres', 0),
                'hydrophobicity': round(float(r.get('hydrophobicity', 0.0)), 3),
                'polar_fraction': round(float(r.get('polar_fraction',  0.0)), 3),
                'box_size':       r.get('box_size', [20.0, 20.0, 20.0]),
                'source':         'polarpocket',
            })
        return result

    # ── Residue profiling ─────────────────────────────────────────────────────
    def _parse_pdb_residues(self) -> list:
        # CRITICAL: use chainA.pdb (author residue numbers) not fixed.pdb
        # (PDBFixer renumbers from 1, breaking ConSurf/B-factor mapping)
        pdb = self._resolve_chainA_pdb() or self._resolve_pdb_path() or self.receptor_pdb
        residues, seen = [], set()
        try:
            with open(pdb) as fh:
                for line in fh:
                    if line[:4] not in ('ATOM', 'HETA'):
                        continue
                    try:
                        res_name = line[17:20].strip()
                        res_num  = int(line[22:26].strip())
                        chain    = line[21]
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        key = (res_name, res_num, chain)
                        if key not in seen:
                            seen.add(key)
                            residues.append({
                                'name': res_name, 'num': res_num,
                                'chain': chain, 'coords': [x, y, z],
                            })
                    except (ValueError, IndexError):
                        continue
        except FileNotFoundError:
            self.log(f'PDB not found for residue parsing: {pdb}', 'WARN')
        return residues

    def _pocket_residue_profile(self, pocket: dict,
                                  radius: float = 6.5) -> dict:
        cx, cy, cz = pocket['center']
        nearby = []
        for res in self._residues:
            rx, ry, rz = res['coords']
            if math.sqrt((rx-cx)**2 + (ry-cy)**2 + (rz-cz)**2) <= radius:
                nearby.append(res['name'])
        n = max(len(nearby), 1)
        # Also collect full residue info (name + number) for annotator
        full_residues = [
            res for res in self._residues
            if math.sqrt(
                (res['coords'][0]-cx)**2 +
                (res['coords'][1]-cy)**2 +
                (res['coords'][2]-cz)**2
            ) <= radius
        ]
        return {
            'n_residues':       n,
            'polar_frac':       round(sum(1 for r in nearby if r in POLAR_AA)    / n, 3),
            'charged_frac':     round(sum(1 for r in nearby if r in CHARGED_AA)  / n, 3),
            'hydrophobic_frac': round(sum(1 for r in nearby if r in HYDRO_AA)    / n, 3),
            'aromatic_frac':    round(sum(1 for r in nearby if r in AROMATIC_AA) / n, 3),
            'metal_coord_n':    sum(1 for r in nearby if r in METAL_AA),
            'top_residues':     [aa for aa, _ in Counter(nearby).most_common(3)],
            'residue_names':    list(set(nearby)),
            'residue_numbers':  list({r['num'] for r in full_residues}),
        }

    # ── Classification ────────────────────────────────────────────────────────
    def _classify(self, pocket: dict, rank: int) -> dict:
        rp   = pocket.get('residue_profile', {})
        vol  = pocket['volume']
        drg  = pocket['druggability']
        dist = pocket.get('dist_from_primary', 0)
        mc   = rp.get('metal_coord_n', 0)
        pf   = rp.get('polar_frac', 0)
        cf   = rp.get('charged_frac', 0)
        hf   = rp.get('hydrophobic_frac', 0)

        if vol < self.MAX_VOL_ION and mc >= self.METAL_CLUSTER_MIN:
            return {'code': 'D', 'label': 'Ion/Metal binding site',
                    'confidence': 'High' if mc >= 3 else 'Moderate'}
        if rank == 0 and drg >= self.MIN_DRUGG_PRIMARY:
            return {'code': 'A', 'label': 'Primary binding / active site',
                    'confidence': 'High'}
        if vol > self.MIN_VOL_PPI and (pf + cf) > 0.5 and drg < 0.55:
            return {'code': 'C', 'label': 'Protein-protein interaction interface',
                    'confidence': 'Moderate'}
        if dist > self.ALLOSTERIC_MIN_DIST and drg >= self.MIN_DRUGG_DRUGGABLE:
            return {'code': 'B', 'label': 'Allosteric / regulatory site',
                    'confidence': 'Moderate'}
        if drg < 0.30 and pf > 0.4 and vol > 200:
            return {'code': 'E', 'label': 'Cryptic / dynamic site',
                    'confidence': 'Low'}
        if hf > 0.6 and drg < 0.30:
            return {'code': 'F', 'label': 'Structural cavity (fold stability)',
                    'confidence': 'Moderate'}
        if drg >= self.MIN_DRUGG_DRUGGABLE:
            return {'code': 'A*', 'label': 'Druggable binding site (unclassified)',
                    'confidence': 'Low'}
        return {'code': 'X', 'label': 'Surface groove (low druggability)',
                'confidence': 'Low'}

    # ── Scientific interpretation ─────────────────────────────────────────────
    def _interpret(self, pocket: dict) -> str:
        code = pocket.get('pocket_type', {}).get('code', 'X')
        rp   = pocket.get('residue_profile', {})
        vol  = pocket['volume']
        drg  = pocket['druggability']
        dist = pocket.get('dist_from_primary', 0)
        top3 = ', '.join(rp.get('top_residues', [])) or 'N/A'

        interp = {
            'A': (
                f"Primary binding/active site (vol={vol:.0f} Å³, drugg={drg:.2f}). "
                f"Dominant residues: {top3}. Highest druggability pocket — canonical "
                f"site for substrate/inhibitor binding. Hydrophobic fraction "
                f"{rp.get('hydrophobic_frac',0):.0%} and aromatic fraction "
                f"{rp.get('aromatic_frac',0):.0%} favour van der Waals and π-stacking. "
                f"Primary candidate for structure-based drug design."
            ),
            'B': (
                f"Allosteric/regulatory site (vol={vol:.0f} Å³, {dist:.1f} Å from primary). "
                f"Dominant residues: {top3}. Distal from the active site — binding here "
                f"can modulate activity through conformational communication without "
                f"competing with the substrate. Polar fraction {rp.get('polar_frac',0):.0%} "
                f"suggests H-bond-mediated allosteric modulation. Selectivity advantage "
                f"over orthosteric inhibition."
            ),
            'C': (
                f"Protein-protein interaction (PPI) interface (vol={vol:.0f} Å³). "
                f"Dominant residues: {top3}. Large shallow surface with charged fraction "
                f"{rp.get('charged_frac',0):.0%} — consistent with electrostatic PPI "
                f"contacts. Disrupting this interface can dissociate pathological "
                f"complexes. Macrocycles or stapled peptides are preferred scaffolds."
            ),
            'D': (
                f"Ion/metal binding site (vol={vol:.0f} Å³, "
                f"{rp.get('metal_coord_n',0)} metal-coord residues: {top3}). "
                f"His/Cys/Asp/Glu cluster indicates coordination of Zn²⁺, Mg²⁺, or Ca²⁺. "
                f"Highly conserved across homologues. Chelating pharmacophores "
                f"(hydroxamate, catechol) can disrupt catalytic metal coordination."
            ),
            'E': (
                f"Cryptic/dynamic site (vol={vol:.0f} Å³, drugg={drg:.2f}). "
                f"Dominant residues: {top3}. Low druggability in static crystal "
                f"structure — pocket transiently accessible through protein dynamics. "
                f"Polar fraction {rp.get('polar_frac',0):.0%} may enable opening upon "
                f"fragment binding. MD simulation recommended to characterise occupancy."
            ),
            'F': (
                f"Structural cavity (vol={vol:.0f} Å³, drugg={drg:.2f}). "
                f"Dominant residues: {top3}. Buried hydrophobic core "
                f"(hydrophobic {rp.get('hydrophobic_frac',0):.0%}) contributes to "
                f"tertiary fold stability. Not a classical drug binding site — "
                f"perturbation risks protein destabilisation."
            ),
        }
        return interp.get(code, interp.get('A*',
            f"Surface pocket (vol={vol:.0f} Å³, drugg={drg:.2f}). "
            f"Dominant residues: {top3}. Druggable but unclassified — "
            f"further characterisation recommended."
        ))

    # ── UniProt annotation ────────────────────────────────────────────────────
    def _fetch_uniprot_annotations(self, pockets: list):
        uid = self.uniprot_id or self._pdb_to_uniprot(self.target_id)
        if not uid:
            return
        self.log(f'Fetching UniProt annotations: {uid}')
        try:
            url  = f'https://www.ebi.ac.uk/proteins/api/features/{uid}'
            resp = requests.get(url, timeout=15,
                                headers={'Accept': 'application/json'})
            if resp.status_code != 200:
                return
            features = resp.json().get('features', [])
            sites    = [f for f in features
                        if f.get('type') in ('ACTIVE_SITE','BINDING','METAL',
                                             'SITE','REGION')]
            self.log(f'  UniProt: {len(sites)} functional annotations')
            for pocket in pockets:
                pocket['uniprot_annotations'] = [
                    {'type': s.get('type',''), 'description': s.get('description','')}
                    for s in sites[:5]
                ]
        except Exception as e:
            self.log(f'UniProt fetch failed: {e}', 'WARN')

    @staticmethod
    def _pdb_to_uniprot(pdb_id: str) -> str:
        try:
            url  = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}'
            data = requests.get(url, timeout=10).json()
            for poly in (data.get('rcsb_entry_container_identifiers', {})
                             .get('polymer_entity_ids', [])):
                ep  = requests.get(
                    f'https://data.rcsb.org/rest/v1/core/polymer_entity'
                    f'/{pdb_id.upper()}/{poly}', timeout=10
                ).json()
                for ref in (ep.get('rcsb_polymer_entity_container_identifiers', {})
                               .get('uniprot_ids') or []):
                    return ref
        except Exception:
            pass
        return ''

    # ── Report ────────────────────────────────────────────────────────────────
    def write_report(self, pockets: list):
        html_out = self.outdir / 'surface_report.html'
        COLORS = {
            'A': '#28a745', 'A*': '#5cb85c', 'B': '#007bff',
            'C': '#fd7e14', 'D': '#6f42c1',  'E': '#ffc107',
            'F': '#6c757d', 'X': '#adb5bd',
        }
        rows = ''
        for p in pockets:
            pt    = p.get('pocket_type', {})
            code  = pt.get('code', 'X')
            color = COLORS.get(code, '#adb5bd')
            rp    = p.get('residue_profile', {})
            ev    = p.get('evidence', {})
            ev_pes   = f"{ev.get('PES', 0):.3f}" if ev else '—'
            ev_label = ev.get('confidence_label', '—') if ev else '—'
            ev_coxtal = (f"✅ {ev['cocrystal_ligand']} ({ev['cocrystal_dist_A']:.1f}Å)"
                         if ev and ev.get('cocrystal_ligand') else '❌')
            ev_cons  = f"{ev.get('consurf_mean','—')}" if ev else '—'
            ev_nma   = f"{ev.get('nma_mobility',0):.3f}" if ev else '—'
            ev_rcn   = f"{ev.get('rcn_centrality',0):.3f}" if ev else '—'
            ev_vina  = f"{ev['best_vina']:.2f}" if ev and ev.get('best_vina') else '—'
            up    = '; '.join(
                f"{a.get('type','')} {a.get('description','')}"
                for a in p.get('uniprot_annotations', [])[:2]
            ) or '—'
            rows += f"""<tr>
  <td><b>{p['rank']}</b></td>
  <td>{p['pocket_id']}</td>
  <td><span style="background:{color};color:white;padding:2px 6px;
      border-radius:4px;font-size:11px">{code}</span>&nbsp;{pt.get('label','')}</td>
  <td>{p['volume']:.0f}</td>
  <td>{p['druggability']:.3f}</td>
  <td>{p.get('druggability_label','')}</td>
  <td>{p.get('dist_from_primary',0):.1f}</td>
  <td>{rp.get('polar_frac',0):.0%}</td>
  <td>{rp.get('hydrophobic_frac',0):.0%}</td>
  <td>{rp.get('charged_frac',0):.0%}</td>
  <td>{rp.get('metal_coord_n',0)}</td>
  <td>{', '.join(rp.get('top_residues',[]))}</td>
  <td>{ev_pes}</td>
  <td>{ev_coxtal}</td>
  <td>{ev_cons}</td>
  <td>{ev_nma}</td>
  <td>{ev_rcn}</td>
  <td>{ev_vina}</td>
  <td style="font-size:11px;max-width:220px">{ev_label}</td>
  <td style="font-size:11px">{up}</td>
</tr>"""

        interp_html = ''
        for p in pockets:
            pt    = p.get('pocket_type', {})
            code  = pt.get('code', 'X')
            color = COLORS.get(code, '#adb5bd')
            interp_html += f"""
<div style="border-left:4px solid {color};padding:10px 16px;margin:10px 0;
            background:#f8f9fa;border-radius:0 8px 8px 0">
  <b>Pocket {p['rank']} (ID {p['pocket_id']})</b>
  <span style="background:{color};color:white;padding:1px 6px;border-radius:4px;
      margin-left:8px;font-size:11px">{pt.get('label','')}</span>
  <small style="color:#666;display:block;margin:4px 0">
    Vol={p['volume']:.0f} Å³ | Drugg={p['druggability']:.3f} |
    {p.get('druggability_label','')} |
    Dist={p.get('dist_from_primary',0):.1f} Å from primary
  </small>
  <p style="margin:6px 0;color:#333;font-size:13px">{p.get('interpretation','')}</p>
</div>"""

        n_druggable = sum(1 for p in pockets
                          if p.get('pocket_type',{}).get('code','X')
                          not in ('F','X','E'))
        html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>SAIDock Surface Analysis — {self.target_id}</title>
<style>
  body{{{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
        padding:20px;max-width:1600px;margin:0 auto}}}}
  h1{{{{color:#212529;border-bottom:3px solid #007bff;padding-bottom:8px}}}}
  h2{{{{color:#495057;margin-top:28px}}}}
  table{{{{border-collapse:collapse;width:100%;font-size:12px}}}}
  th{{{{background:#343a40;color:#fff;padding:8px 6px;position:sticky;top:0}}}}
  td{{{{padding:6px;border-bottom:1px solid #dee2e6;vertical-align:top}}}}
  tr:hover td{{{{background:#f1f3f5}}}}
  .badge{{{{display:inline-block;padding:2px 8px;border-radius:4px;
          color:white;font-size:11px;margin:2px}}}}
</style>
  </head><body>
<h1>🔬 SAIDock Surface Analysis — {self.target_id}</h1>
<p>
  <b>{len(pockets)}</b> pockets analysed &nbsp;|&nbsp;
  <b>{n_druggable}</b> biologically relevant (types A/B/C/D) &nbsp;|&nbsp;
  Source: {pockets[0].get('source','N/A') if pockets else 'N/A'}
</p>
<div>
  <span class="badge" style="background:#28a745">A: Primary site</span>
  <span class="badge" style="background:#007bff">B: Allosteric</span>
  <span class="badge" style="background:#fd7e14">C: PPI interface</span>
  <span class="badge" style="background:#6f42c1">D: Ion/Metal</span>
  <span class="badge" style="background:#ffc107;color:#333">E: Cryptic</span>
  <span class="badge" style="background:#6c757d">F: Structural</span>
</div>
<h2>Pocket Atlas</h2>
<div style="overflow-x:auto">
<table>
<thead><tr>
  <th>#</th><th>ID</th><th>Type</th><th>Vol (Å³)</th><th>Drugg.</th>
  <th>Assessment</th><th>Dist (Å)</th><th>Polar%</th><th>Hydro%</th>
  <th>Charged%</th><th>Metal AA</th><th>Top Residues</th>
  <th>PES</th><th>CoXtal</th><th>ConSurf</th><th>NMA mob.</th>
  <th>RCN cent.</th><th>Best Vina</th><th>Evidence</th><th>UniProt</th>
</tr></thead>
<tbody>{rows}</tbody>
</table></div>
<h2>Scientific Interpretation</h2>
{interp_html}
<hr><small style="color:#999">SAIDock v1.0.0 SurfaceAnalyser</small>
</body></html>"""

        # ── NGL 3D viewer post-processing (no f-string) ──────────────
        try:
            import json as _j
            _pd = []
            for _p in pockets:
                _c = (_p.get('center') or _p.get('centroid') or
                      _p.get('coords') or [0, 0, 0])
                if len(_c) >= 3:
                    _ev = _p.get('evidence', {})
                    _pd.append({
                        'id':   str(_p.get('pocket_id', _p.get('id', '?'))),
                        'cx':   float(_c[0]),
                        'cy':   float(_c[1]),
                        'cz':   float(_c[2]),
                        'dtss': round(float(
                            _ev.get('PES', _p.get('PES', 0))), 3),
                        'conf': str(_ev.get(
                            'confidence_label',
                            _p.get('confidence_label', ''))),
                    })
            _rp = str(self.outdir /
                      f'{self.target_id}_chainA.pdb').replace("\\", "/")
            _ngl_block = (_NGL_BODY_TPL
                         .replace("__POCKETJSON__", _j.dumps(_pd))
                         .replace("__RECPATH__",    _rp))
            html = html.replace(
                '</head>', _NGL_HEAD_SCRIPT + '</head>', 1)
            html = html.replace('</body>', _ngl_block + '</body>', 1)
        except Exception as _ngl_err:
            import traceback as _tb
            (self.outdir / "ngl_error.log").write_text(_tb.format_exc())
        # ─────────────────────────────────────────────────────────────
        html_out.write_text(html)
        self.log(f'Report: {html_out}', 'OK')

    # ── Utilities ─────────────────────────────────────────────────────────────
    @staticmethod
    def _dist(a, b) -> float:
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))

    @staticmethod
    def _drugg_label(score: float) -> str:
        if score >= 0.7: return '🟢 Highly druggable'
        if score >= 0.5: return '🟡 Druggable'
        if score >= 0.3: return '🟠 Moderately druggable'
        return '🔴 Poorly druggable'

    @staticmethod
    def _write_csv(pockets: list, path: Path):
        import csv
        if not pockets:
            return
        rows = []
        for p in pockets:
            rp = p.get('residue_profile', {})
            pt = p.get('pocket_type', {})
            rows.append({
                'rank':              p['rank'],
                'pocket_id':         p['pocket_id'],
                'type_code':         pt.get('code', ''),
                'type_label':        pt.get('label', ''),
                'confidence':        pt.get('confidence', ''),
                'center_x':          round(p['center'][0], 3),
                'center_y':          round(p['center'][1], 3),
                'center_z':          round(p['center'][2], 3),
                'volume_A3':         p['volume'],
                'druggability':      p['druggability'],
                'druggability_label':p.get('druggability_label', ''),
                'dist_from_primary': round(p.get('dist_from_primary', 0), 2),
                'n_residues':        rp.get('n_residues', 0),
                'polar_frac':        rp.get('polar_frac', 0),
                'hydrophobic_frac':  rp.get('hydrophobic_frac', 0),
                'charged_frac':      rp.get('charged_frac', 0),
                'aromatic_frac':     rp.get('aromatic_frac', 0),
                'metal_coord_n':     rp.get('metal_coord_n', 0),
                'top_residues':      ','.join(rp.get('top_residues', [])),
                'n_spheres':         p.get('n_spheres', 0),
                'hydrophobicity':    p.get('hydrophobicity', 0),
                'source':            p.get('source', ''),
                'interpretation':    p.get('interpretation', '')[:400],
            })
        with open(path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
