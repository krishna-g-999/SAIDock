
def _is_float(s):
    try: float(s); return True
    except (ValueError, TypeError): return False

#!/usr/bin/env python3
"""
SAIDock PocketAnnotator v2 — Multi-evidence pocket role determination.

Evidence levels:
  L1  Structural : co-crystal HETATM + other-chain peptide/ligand detection
  L2  Database   : PDBe Binding Sites API (EBI) — verified binding residues
  L3  Evolution  : ConSurfDB conservation scores (3 URL patterns tried)
  L4a Dynamic    : Numpy ANM / ProDy normal mode fluctuations
  L4b Network    : Residue Contact Network betweenness centrality
  L5  Docking    : Auto-load SAIDock batch results from results/ directory

Dynamic scoring — scientifically correct:
  PRIMARY site   = rigid (low NMA) + well-connected (high RCN) = active site hub
  ALLOSTERIC     = mobile (high NMA) + well-connected = conformational hinge
  dynamic_score  = max(primary_dynamic, allosteric_dynamic)
  pocket_role_hint stored per pocket to guide classification
"""
import shutil
import csv, json, math, traceback, warnings
from pathlib import Path
from typing   import Dict, List, Optional
import numpy as np
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

def _make_session(retries: int = 2, timeout: int = 5) -> requests.Session:
    """Session with retry + short timeout — fail fast on dead servers."""
    session = requests.Session()
    adapter = HTTPAdapter(
        max_retries=Retry(total=retries, backoff_factor=0.3,
                          status_forcelist=[429, 500, 502, 503])
    )
    session.mount('https://', adapter)
    session.mount('http://',  adapter)
    return session

_SESSION = _make_session()

warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=UserWarning)

METALS   = {'ZN','MG','CA','FE','MN','CO','NI','CU','NA','K',
            'FE2','FE3','MO','W','V','SE'}
SKIP_HET = {'HOH','WAT','DOD','SO4','PO4','GOL','EDO','PEG','MPD',
            'BME','DMS','ACT','FMT','ACE','MES','HEPES','TRS',
            'IOD','CL','BR','IOD','SO3','NO3','NO2','PGE',
            'BU3','IPA','ETH','DIO','THF','ACN'}

# Biologically relevant ligands — co-crystal with these = CONFIRMED binding site
# regardless of PES threshold. Add new entries as needed.
BIOLOGICAL_LIGANDS = {
    # Peptide/peptidomimetic ligands
    '41P','ETG','NRF','NQO','KBD',
    # Kinase ATP-site ligands
    'ATP','ADP','AMP','ANP','ACP','STU','SB1','PD1','CDK',
    # Common pharmacophore reference ligands
    'NAD','FAD','FMN','HEM','HEC','MG3','SF4','FES',
    # NF-kB / Nrf2 pathway
    'SFN','CUR','RES','QUE','KAE',
}

# Evidence weights — must sum to 1.0
WEIGHTS = dict(structural=0.30, conservation=0.25,
               dynamic=0.25, database=0.10, docking=0.10)



def _pocket_conservation(cons_scores: dict, pocket_residues: list,
                          fallback_mean: float = 0.5) -> float:
    """Average conservation score over residues lining this specific pocket.

    If pocket_residues is empty (L2 returned nothing), falls back to the
    mean over the nearest 10 residues by index — still pocket-specific.
    """
    if not cons_scores:
        return fallback_mean

    if pocket_residues:
        vals = [cons_scores[r] for r in pocket_residues if r in cons_scores]
        if vals:
            return round(sum(vals) / len(vals), 4)

    # Fallback: use global mean (only if no residue info available)
    all_vals = list(cons_scores.values())
    return round(sum(all_vals) / len(all_vals), 4) if all_vals else fallback_mean


def _pocket_cons_mean(cons_scores: dict, pocket_residues: list) -> float:
    """Return conservation mean over pocket-lining residues only.
    
    If pocket_residues is populated (from L2 BioLiP or HETATM contacts),
    averages only those residues. Otherwise falls back to protein mean.
    For ConSurf grades (1–9 scale) returns raw mean; caller normalizes.
    For HMMER (0–1 scale) returns normalized mean directly.
    """
    if not cons_scores:
        # Return ConSurf-scale 5.0 (midpoint) when no data
        return 5.0 if max(cons_scores.values(), default=0) > 1.0 else 0.5

    # Pocket-specific: prefer residues known to contact the ligand
    if pocket_residues:
        vals = [cons_scores[r] for r in pocket_residues if r in cons_scores]
        if len(vals) >= 3:   # need at least 3 residues for meaningful mean
            return sum(vals) / len(vals)

    # Fallback: global protein mean (less informative but never crashes)
    all_vals = list(cons_scores.values())
    return sum(all_vals) / len(all_vals) if all_vals else 0.5

class PocketAnnotator:
    """Adds multi-evidence annotations to pocket dicts from SurfaceAnalyser."""

    def __init__(self, pdb_id: str, raw_pdb_path: str,
                 fixed_pdb_path: str, chain: str = 'A',
                 logger=None, cache_dir: str = '.saidock_cache'):
        self.pdb_id     = pdb_id.upper()
        self.raw_pdb    = Path(raw_pdb_path)
        self.fixed_pdb  = Path(fixed_pdb_path)
        self.chain      = chain
        self.logger     = logger
        self.cache_dir  = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._hetatm_entries : List[Dict] = []    # {name,x,y,z,is_metal,is_organic}
        self._other_chains   : List[Dict] = []    # ATOM records from non-target chains
        self._pdbe_res       : set        = set() # PDBe binding residue numbers
        self._consurf        : Dict[int, float] = {}
        self._nma_fluct      : Dict[int, float] = {}
        self._rcn_cent       : Dict[int, float] = {}
        self._docking_results: List[Dict] = []

    def _get_chainA_pdb(self) -> Path:
        """
        Return the chainA PDB (author residue numbers preserved).
        This is the SINGLE SOURCE OF TRUTH for all residue numbering.
        fixed.pdb renumbers from 1 — never use it for residue lookups.
        Priority: *_chainA.pdb → raw_pdb → fixed_pdb (last resort only)
        """
        for parent in [self.raw_pdb.parent,
                       self.fixed_pdb.parent,
                       self.fixed_pdb.parent.parent]:
            hits = sorted(parent.glob('*_chainA.pdb'))
            if hits:
                return hits[0]
        return self.raw_pdb if self.raw_pdb.exists() else self.fixed_pdb

    def log(self, msg: str, level: str = 'INFO'):
        if self.logger:
            self.logger.log(msg, level)
        else:
            print(f'  [{level}] {msg}')

    # ── Main entry ─────────────────────────────────────────────────────────
    def annotate_all(self, pockets: list, docking_results: list = None) -> list:
        if docking_results:
            self._docking_results = docking_results
        else:
            self._auto_load_docking()   # L5: scan results/ directory

        self.log('Loading evidence layers...')
        self._scan_hetatm_and_chains()  # L1
        self._query_pdbe_binding()      # L2
        self._load_consurf()            # L3
        self._compute_nma()             # L4a
        self._compute_rcn()             # L4b

        for pocket in pockets:
            self._annotate_pocket(pocket)
        self._reclassify(pockets)
        pockets = self._merge_subpockets(pockets)
        self.log(f'Evidence annotation complete — {len(pockets)} pockets', 'OK')
        return pockets

    # ── L1: HETATM + other chains ───────────────────────────────────────────
    def _scan_hetatm_and_chains(self):
        """
        Scan raw PDB for:
          a) HETATM records (small molecule ligands + metals)
          b) ATOM records in OTHER chains (co-crystal peptides, protein partners)
        Both are biological evidence for nearby pocket function.
        """
        pdb = self.raw_pdb if self.raw_pdb.exists() else self.fixed_pdb
        if not pdb.exists():
            self.log(f'HETATM scan: PDB not found', 'WARN')
            return

        organic_n, metal_n, other_chain_n = 0, 0, 0
        with open(pdb) as fh:
            for line in fh:
                rec = line[:6].strip()
                # HETATM scan
                if rec == 'HETATM':
                    try:
                        name   = line[17:20].strip()
                        alt    = line[16].strip()
                        if alt and alt not in ('A', ' ', ''):
                            continue
                        if name in SKIP_HET:
                            continue
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        is_metal  = name.upper() in METALS
                        is_organic = not is_metal
                        self._hetatm_entries.append({
                            'name': name, 'x': x, 'y': y, 'z': z,
                            'is_metal': is_metal, 'is_organic': is_organic,
                        })
                        if is_metal: metal_n += 1
                        else:        organic_n += 1
                    except (ValueError, IndexError):
                        continue

                # Other-chain ATOM scan (co-crystal peptides, protein partners)
                elif rec == 'ATOM':
                    try:
                        chain_id = line[21]
                        if chain_id == self.chain:
                            continue
                        atom_name = line[12:16].strip()
                        if atom_name != 'CA':   # Cα only for speed
                            continue
                        res_name = line[17:20].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        self._other_chains.append({
                            'name': res_name, 'chain': chain_id,
                            'x': x, 'y': y, 'z': z,
                            'is_metal': False, 'is_organic': True,
                        })
                        other_chain_n += 1
                    except (ValueError, IndexError):
                        continue

        self.log(f'L1: {organic_n} organic ligands + {metal_n} metals + '
                 f'{other_chain_n} other-chain residues')

    def _closest_ligand(self, center: list, radius: float = 8.0) -> Dict:
        """Return closest HETATM or other-chain residue within radius."""
        cx, cy, cz = center
        best_d, best = radius + 1, None
        all_hits = self._hetatm_entries + self._other_chains
        for h in all_hits:
            d = math.sqrt((h['x']-cx)**2 + (h['y']-cy)**2 + (h['z']-cz)**2)
            if d < best_d:
                best_d, best = d, h
        if best and best_d <= radius:
            return {
                'name':       best['name'],
                'dist':       round(best_d, 2),
                'is_metal':   best['is_metal'],
                'source':     'hetatm' if best in self._hetatm_entries
                              else 'other_chain',
            }
        return {}

    # ── L2: PDBe Binding Sites API ──────────────────────────────────────────
    def _query_pdbe_binding(self):
        """
        PDBe API: https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{PDB_ID}
        Returns experimentally identified binding sites + residue lists.
        Falls back to RCSB struct_site if PDBe fails.
        """
        cache_file = self.cache_dir / f'pdbe_binding_{self.pdb_id}.json'
        if cache_file.exists():
            data = json.loads(cache_file.read_text())
            self._pdbe_res = set(data.get('residue_numbers', []))
            self.log(f'L2 PDBe (cached): {len(self._pdbe_res)} binding residues')
            return

        res_nums = set()

        # Try PDBe first
        try:
            url  = (f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/'
                    f'{self.pdb_id.lower()}')
            resp = _SESSION.get(url, timeout=5,
                                headers={'Accept': 'application/json'})
            if resp.status_code == 200:
                data = resp.json()
                # PDBe graph-api returns list of dicts per ligand
                raw_data = data if isinstance(data, list) else                            data.get(self.pdb_id.lower(),
                           data.get(self.pdb_id.upper(), []))
                if isinstance(raw_data, dict):
                    raw_data = [raw_data]
                for entry in raw_data:
                    # Handle both binding_sites and bound_ligand_residues formats
                    for key in ('site_residues','residues','binding_residues',
                                'contact_residues'):
                        for res in entry.get(key, []):
                            ch = (res.get('chain_id') or
                                  res.get('struct_asym_id') or
                                  res.get('chain') or self.chain)
                            if ch != self.chain:
                                continue
                            rn = (res.get('residue_number') or
                                  res.get('author_residue_number') or
                                  res.get('seq_id') or
                                  res.get('label_seq_id'))
                            if rn:
                                try: res_nums.add(int(rn))
                                except (ValueError, TypeError): pass
                if res_nums:
                    self.log(f'L2 PDBe: {len(res_nums)} binding residues', 'OK')
                else:
                    self.log('L2 PDBe: 0 residues found in response', 'WARN')
        except (requests.exceptions.Timeout,
                requests.exceptions.ConnectionError) as e:
            self.log(f'PDBe timeout/connection: {e}', 'WARN')
        except Exception as e:
            self.log(f'PDBe binding sites failed: {e}', 'WARN')

        # RCSB struct_site fallback
        if not res_nums:
            try:
                url  = (f'https://data.rcsb.org/rest/v1/core/'
                        f'entry/{self.pdb_id.upper()}')
                resp = requests.get(url, timeout=5)  # no retry — fail fast on dead ConSurf server
                if resp.status_code == 200:
                    data = resp.json()
                    for site in data.get('struct_site', []):
                        for detail in site.get('rcsb_struct_site_gen', []):
                            if detail.get('label_asym_id') == self.chain:
                                rn = detail.get('label_seq_id') or \
                                     detail.get('auth_seq_id')
                                if rn:
                                    try: res_nums.add(int(rn))
                                    except (ValueError, TypeError): pass
                    self.log(f'L2 RCSB fallback: {len(res_nums)} binding residues')
            except Exception as e:
                self.log(f'RCSB binding sites fallback failed: {e}', 'WARN')

        self._pdbe_res = res_nums
        if res_nums:  # only cache non-empty results
                cache_file.write_text(json.dumps(
            {'residue_numbers': list(res_nums)}, indent=2))

    # ── L3: ConSurfDB conservation ───────────────────────────────────────────
    def _load_consurf(self):
        """
        ConSurfDB pre-computed conservation scores (1=variable, 9=conserved).
        Tries 3 known URL patterns. Falls back to SIFTS + UniRef alignment.
        """
        cache_file = self.cache_dir / f'consurf_{self.pdb_id}_{self.chain}.json'
        if cache_file.exists():
            self._consurf = {int(k): v
                             for k, v in json.loads(cache_file.read_text()).items()}
            self.log(f'L3 ConSurf (cached): {len(self._consurf)} residue scores')
            return

        pid = self.pdb_id.lower()
        CH  = self.chain

        url_patterns = [
            # Pattern 1: most common ConSurfDB format (4-char PDB, 1-char chain)
            f'https://consurf.tau.ac.il/ConSurfDB/data/{pid}/{CH}/consurf_grades.txt',
            # Pattern 2: subdirectory by middle chars
            f'https://consurf.tau.ac.il/ConSurfDB/data/{pid[1:3]}/{pid}/{CH}/consurf_grades.txt',
            # Pattern 3: flat layout
            f'https://consurf.tau.ac.il/ConSurfDB/{pid}/{CH}/consurf_grades.txt',
            # Pattern 4: precomputed grades
            f'https://consurf.tau.ac.il/results/precomputed/{pid.upper()}_{CH}.grades',
        ]

        # ── Priority 0: HMMER + Shannon entropy ──────────────────────
        scores = self._run_hmmer_conservation()
        if scores:
            self._consurf = scores   # ← store, not return
            return
        # ── Priority 1: local ConSurf grades file ─────────────────────
        for url in url_patterns:
            try:
                resp = requests.get(url, timeout=5)  # 8s max — fail fast
                if resp.status_code == 200 and 'GRADE' in resp.text.upper():
                    scores = self._parse_consurf_grades(resp.text)
                    if len(scores) >= 10:
                        self._consurf = scores
                        cache_file.write_text(json.dumps(scores, indent=2))
                        self.log(f'L3 ConSurf: {len(scores)} residues', 'OK')
                        return
            except Exception:
                continue

        self.log('L3 ConSurfDB not reachable — using B-factor gradient proxy', 'WARN')
        self._bfactor_gradient_proxy()

    def _parse_consurf_grades(self, text: str) -> Dict[int, float]:
        """Parse ConSurf grades file (multiple version formats supported)."""
        scores = {}
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('POS'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                pos = int(parts[0])
                # Score is usually in column 3 or 4; look for a 1-9 float
                for p in parts[2:6]:
                    p2 = p.strip('*[]')
                    try:
                        val = float(p2)
                        if 1.0 <= val <= 9.0:
                            scores[pos] = val
                            break
                    except ValueError:
                        continue
            except (ValueError, IndexError):
                continue
        return scores

    def _run_hmmer_conservation(self) -> dict:
        """L3 Priority 1: phmmer + Shannon entropy → conservation per PDB residue.

        Returns dict keyed by PDB residue number (matching fpocket _lining_resnums),
        NOT by sequential HMMER position. This is the critical remapping step.
        Reference: Valdar (2002) Proteins. PMID:11835497
        """
        import subprocess, tempfile as _tf
        if not shutil.which('phmmer'):
            return {}
        db_candidates = [
            Path('data/seqdb/pdb_seqres.fasta'),
            Path('/mirror/pdb_seqres.fasta'),
            Path('/data/pdb_seqres.fasta'),
        ]
        db = next((p for p in db_candidates
                   if p.exists() and p.stat().st_size > 100_000), None)
        if not db:
            self.log('L3 HMMER: pdb_seqres.fasta not found', 'WARN')
            return {}

        fasta, pdb_resnums = self._extract_fasta_seq()
        if not fasta or not pdb_resnums:
            return {}

        try:
            with _tf.TemporaryDirectory() as tmp:
                tmp     = Path(tmp)
                fa      = tmp / f'{self.pdb_id}.fasta'
                sto     = tmp / f'{self.pdb_id}.sto'
                fa.write_text(fasta)
                ret = subprocess.run(
                    ['phmmer', '-A', str(sto), '--incE', '0.001',
                     '--cpu', '4', str(fa), str(db)],
                    capture_output=True, text=True, timeout=300
                )
                if not sto.exists() or sto.stat().st_size < 100:
                    self.log(f'L3 phmmer: no hits (exit {ret.returncode})', 'WARN')
                    return {}

                # scores keyed by sequential position 1..len(query)
                pos_scores = self._shannon_from_sto(sto)
                if not pos_scores:
                    return {}

                # ── CRITICAL REMAP: sequential pos → PDB residue number ──────
                # pos_scores[1] = conservation of 1st residue in FASTA
                # pdb_resnums[0] = PDB residue number of 1st residue
                # These can differ if PDB starts at residue 5, has gaps, etc.
                resnum_scores = {
                    pdb_resnums[pos - 1]: score
                    for pos, score in pos_scores.items()
                    if 0 <= (pos - 1) < len(pdb_resnums)
                }

                self.log(
                    f'L3 HMMER+Shannon: {len(resnum_scores)} residues '
                    f'(PDB resnums {min(resnum_scores)}-{max(resnum_scores)}) '
                    f'| db={db.name}', 'INFO'
                )
                return resnum_scores

        except Exception as e:
            self.log(f'L3 HMMER error: {e}', 'WARN')
        return {}


    def _shannon_from_sto(self, sto_file: Path) -> dict:
        """Shannon entropy from Stockholm MSA → per-QUERY-RESIDUE conservation.

        Key fix: MSA has N columns (N ≥ query length due to inserted gaps from
        homologs). We parse the query row and only score non-gap query positions,
        mapping column index → query residue number correctly.
        """
        import math
        seqs    = {}
        order   = []
        with open(sto_file) as fh:
            for line in fh:
                s = line.rstrip()
                if not s or s.startswith('#') or s == '//':
                    continue
                parts = s.split()
                if len(parts) == 2:
                    name, seq = parts
                    if name not in seqs:
                        order.append(name)
                    seqs[name] = seqs.get(name, '') + seq

        if not seqs or not order:
            return {}

        # Query is always the FIRST sequence in phmmer Stockholm output
        query_seq = seqs[order[0]]
        seq_list  = list(seqs.values())
        AA        = set('ACDEFGHIKLMNPQRSTVWY')
        H_max     = math.log2(20)
        scores    = {}
        query_pos = 0  # tracks position in the UNGAPPED query sequence

        for col in range(len(query_seq)):
            q_aa = query_seq[col].upper()
            is_gap = q_aa in ('.', '-', '*', 'X', ' ')

            if not is_gap:
                query_pos += 1   # this column maps to residue query_pos
                counts = {}
                for seq in seq_list:
                    if col < len(seq):
                        aa = seq[col].upper()
                        if aa in AA:
                            counts[aa] = counts.get(aa, 0) + 1
                total = sum(counts.values())
                if total > 0:
                    H = -sum((c/total) * math.log2(c/total)
                             for c in counts.values() if c > 0)
                    scores[query_pos] = round(1.0 - H / H_max, 4)

        return scores


    def _extract_fasta_seq(self):
        """Extract Cα sequence from chainA PDB (author residue numbers).

        CRITICAL: Must use chainA.pdb (author numbering) NOT fixed.pdb
        (which MODELLER renumbers 1..N). This ensures HMMER pos→resnum
        mapping matches _lining_resnums exactly.

        Returns (fasta_str, ordered_author_resnums)
        """
        AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
               'GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
               'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
               'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
        seq, seen, resnums = [], set(), []

        # ── MUST use chainA.pdb (author numbering) ──────────────────────
        # fixed.pdb is renumbered 1..N by MODELLER → never matches
        # _lining_resnums which come from chainA author numbering
        pdb_file = None
        try:
            chain_pdb = self._get_chainA_pdb()
            if chain_pdb and Path(str(chain_pdb)).exists():
                pdb_file = Path(str(chain_pdb))
        except Exception:
            pass

        # Fallback chain: any *_chainA.pdb file
        if pdb_file is None:
            hits = (list(self.results_dir.glob(f'{self.pdb_id}_chainA.pdb')) +
                    list(self.results_dir.glob(f'{self.pdb_id.lower()}_chainA.pdb')))
            if hits:
                pdb_file = hits[0]

        # Last resort: fixed_pdb (numbering will be wrong, but better than nothing)
        if pdb_file is None:
            for attr in ('fixed_pdb', 'receptor_pdb', 'prepared_pdb'):
                c = getattr(self, attr, None)
                if c and Path(str(c)).exists():
                    pdb_file = Path(str(c))
                    self.log(
                        f'L3 WARNING: using {Path(str(c)).name} for FASTA — '
                        f'residue numbering may not match _lining_resnums', 'WARN'
                    )
                    break

        if pdb_file is None:
            return '', []

        try:
            chain = getattr(self, 'chain', 'A')
            for line in pdb_file.read_text().splitlines():
                if not line.startswith('ATOM'):
                    continue
                if line[12:16].strip() != 'CA':
                    continue
                # Only target chain
                if line[21:22].strip() not in ('', chain):
                    continue
                try:
                    rnum  = int(line[22:26].strip())
                    rname = line[17:20].strip()
                except (ValueError, IndexError):
                    continue
                if rnum not in seen:
                    seen.add(rnum)
                    seq.append(AA3.get(rname, 'X'))
                    resnums.append(rnum)

            if seq:
                self.log(
                    f'L3 FASTA from {pdb_file.name}: {len(seq)} residues '
                    f'(resnums {min(resnums)}–{max(resnums)})'
                )
                return f'>{self.pdb_id}\n{"".join(seq)}\n', resnums
        except Exception as e:
            self.log(f'L3 _extract_fasta_seq error: {e}', 'WARN')

        return '', []


    def _bfactor_gradient_proxy(self):
        """
        B-factor proxy with SPATIAL GRADIENT normalisation:
        Instead of global normalisation (which makes all scores ~9 for
        well-diffracting crystals), use local percentile-based scoring.
        Residues in the BOTTOM 20% of B-factors = most rigid = score 9.
        Residues in TOP 20% = most flexible = score 1.
        """
        # CRITICAL: use raw/chainA PDB — PDBFixer sets all B-factors to
        # uniform values, collapsing variance and making proxy useless.
        # Priority: raw_pdb (original download) > chainA > fixed
        pdb_candidates = [
            self.raw_pdb,
            self.raw_pdb.parent / (self.raw_pdb.stem + '_chainA.pdb'),
            self.fixed_pdb,
        ]
        pdb = next((p for p in pdb_candidates if p.exists()), None)
        if not pdb:
            return
        bfactors = {}
        with open(pdb) as fh:
            for line in fh:
                if line[:4] != 'ATOM' or line[12:16].strip() != 'CA':
                    continue
                try:
                    chain_col = line[21]
                    if chain_col != self.chain:
                        continue
                    rn   = int(line[22:26])
                    bfac = float(line[60:66])
                    bfactors[rn] = bfac
                except (ValueError, IndexError):
                    continue
        if not bfactors:
            return
        vals     = np.array(list(bfactors.values()))
        # Rank-based: avoid global min-max collapse for uniform B-factors
        from scipy.stats import rankdata
        try:
            ranks    = rankdata(vals)           # ascending: low B = low rank
            n        = len(vals)
            # Low B (rigid) → high conservation score (9); High B → score (1)
            norm     = 1.0 - (ranks - 1) / (n - 1 + 1e-9)  # 1=most rigid
            score_arr = 1.0 + 8.0 * norm       # maps to 1-9
            for (rn, _), sc in zip(sorted(bfactors.items()), score_arr):
                self._consurf[rn] = round(float(sc), 2)
            self.log(f'L3 B-factor gradient proxy: {len(self._consurf)} residues '
                     f'(range {score_arr.min():.1f}–{score_arr.max():.1f})')
        except ImportError:
            # scipy not available — simple percentile
            p20 = np.percentile(vals, 20)
            p80 = np.percentile(vals, 80)
            rng = max(p80 - p20, 1e-6)
            for rn, bf in bfactors.items():
                norm = max(0.0, min(1.0, (bf - p20) / rng))
                self._consurf[rn] = round(9.0 - 8.0 * norm, 2)

    # ── L4a: NMA (numpy ANM, ProDy optional) ───────────────────────────────
    def _compute_nma(self):
        """Anisotropic Network Model: per-residue mean-square fluctuation."""
        # MUST use chainA.pdb — same numbering as lining_resnums
        pdb = self._get_chainA_pdb()
        if not pdb or not pdb.exists():
            self.log('NMA: chainA PDB not found', 'WARN')
            return

        ca_coords, ca_resnums = [], []
        with open(pdb) as fh:
            for line in fh:
                if line[:4] != 'ATOM' or line[12:16].strip() != 'CA':
                    continue
                try:
                    ch = line[21]
                    if ch != self.chain:
                        continue
                    ca_resnums.append(int(line[22:26]))
                    ca_coords.append([float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])])
                except (ValueError, IndexError):
                    continue

        if len(ca_coords) < 5:
            return

        # Try ProDy (faster for large structures)
        try:
            import prody as _prody
            _prody.confProDy(verbosity='none')
            protein = _prody.parsePDB(str(pdb), chain=self.chain, subset='calpha')
            if protein is not None and len(protein) >= 5:
                anm = _prody.ANM(self.pdb_id)
                anm.buildHessian(protein, cutoff=13.0)
                anm.calcModes(n_modes=min(20, len(protein)*3-7))
                sqf      = _prody.calcSqFlucts(anm)
                resnums  = protein.getResnums()
                sqf_norm = (sqf - sqf.min()) / (sqf.max() - sqf.min() + 1e-9)
                self._nma_fluct = {int(rn): round(float(sf), 4)
                                   for rn, sf in zip(resnums, sqf_norm)}
                self.log(f'L4a NMA (ProDy): {len(self._nma_fluct)} residues', 'OK')
                return
        except Exception:
            pass

        # Numpy ANM fallback
        self.log('L4a NMA (numpy ANM)')
        coords  = np.array(ca_coords, dtype=np.float64)
        n       = len(coords)
        cutoff  = 13.0
        N3      = n * 3
        H = np.zeros((N3, N3), dtype=np.float64)
        for i in range(n):
            for j in range(i + 1, n):
                diff  = coords[j] - coords[i]
                dist2 = np.dot(diff, diff)
                if dist2 > cutoff ** 2:
                    continue
                block = -(1.0 / dist2) * np.outer(diff, diff)
                ii, jj = i*3, j*3
                H[ii:ii+3, jj:jj+3]  = block
                H[jj:jj+3, ii:ii+3]  = block
                H[ii:ii+3, ii:ii+3] -= block
                H[jj:jj+3, jj:jj+3] -= block
        try:
            eigvals, eigvecs = np.linalg.eigh(H)
        except np.linalg.LinAlgError:
            self.log('NMA: eigen-decomposition failed', 'WARN')
            return
        nz  = np.where(eigvals > 1e-6)[0]
        sqf = np.zeros(n)
        for k in nz[:20]:
            v = eigvecs[:, k]
            for i in range(n):
                sqf[i] += (v[i*3]**2 + v[i*3+1]**2 + v[i*3+2]**2) / eigvals[k]
        sqf_min, sqf_max = sqf.min(), sqf.max() + 1e-9
        self._nma_fluct = {int(rn): round(float((s-sqf_min)/(sqf_max-sqf_min)), 4)
                           for rn, s in zip(ca_resnums, sqf)}
        self.log(f'L4a NMA (numpy): {n} residues', 'OK')

    # ── L4b: Residue Contact Network ────────────────────────────────────────
    def _compute_rcn(self):
        """Cα contact network: betweenness centrality = allosteric hub score."""
        try:
            import networkx as nx
        except ImportError:
            self.log('networkx not available — skip RCN', 'WARN')
            return

        # MUST use chainA.pdb — same numbering as lining_resnums
        pdb = self._get_chainA_pdb()
        if not pdb or not pdb.exists():
            self.log('RCN: chainA PDB not found', 'WARN')
            return

        ca_coords, ca_resnums = [], []
        with open(pdb) as fh:
            for line in fh:
                if line[:4] != 'ATOM' or line[12:16].strip() != 'CA':
                    continue
                try:
                    ch = line[21]
                    if ch != self.chain:
                        continue
                    ca_resnums.append(int(line[22:26]))
                    ca_coords.append([float(line[30:38]),
                                      float(line[38:46]),
                                      float(line[46:54])])
                except (ValueError, IndexError):
                    continue

        if len(ca_coords) < 5:
            return

        coords = np.array(ca_coords)
        n      = len(coords)
        G      = nx.Graph()
        G.add_nodes_from(ca_resnums)
        for i in range(n):
            for j in range(i+1, n):
                d = np.linalg.norm(coords[i] - coords[j])
                if d < 7.0:
                    G.add_edge(ca_resnums[i], ca_resnums[j],
                               weight=1.0/(d+1e-6))

        cent  = nx.betweenness_centrality(G, normalized=True)
        c_arr = np.array(list(cent.values()))
        c_min, c_max = c_arr.min(), c_arr.max() + 1e-9
        self._rcn_cent = {k: round(float((v-c_min)/(c_max-c_min)), 4)
                          for k, v in cent.items()}
        self.log(f'L4b RCN: {G.number_of_nodes()} nodes, '
                 f'{G.number_of_edges()} edges', 'OK')

    # ── L5: Auto-load docking results ───────────────────────────────────────
    def _auto_load_docking(self):
        """
        Scan results/ directory for batch CSVs containing docking results
        for this target (matched by PDB ID in filename or Target column).
        """
        results_dir = Path('results')
        if not results_dir.exists():
            return
        found = []
        for csv_path in results_dir.rglob('batch_results.csv'):
            try:
                with open(csv_path) as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        target = (row.get('Target') or row.get('target') or
                                  row.get('pdb_id') or '')
                        if self.pdb_id.upper() not in target.upper():
                            continue
                        # Try to find pocket center from docking output dir
                        # Map all known SAIDock CSV column names
                        _SCORE_COLS = ['best_score','Vina_best','vina_best',
                                       'Best_Vina','best_vina','vina_score',
                                       'Score','score']
                        score_raw = next((row[c] for c in _SCORE_COLS
                                          if c in row and row[c]), 0)
                        try:
                            score = float(score_raw) if score_raw else 0.0
                        except (ValueError, TypeError):
                            score = 0.0
                        found.append({'score': score, 'row': row,
                                      'pocket_center': None,
                                      'is_global': True})  # no pocket-level coords
            except Exception:
                continue

        if found:
            self._docking_results = found
            self.log(f'L5 Docking: loaded {len(found)} results for {self.pdb_id}')

    # ── Per-pocket annotation ───────────────────────────────────────────────
    def _annotate_pocket(self, pocket: dict):
        center   = pocket['center']
        res_nums = pocket.get('_lining_resnums', [])

        # L1: structural evidence
        lig_hit    = self._closest_ligand(center, radius=6.0)
        has_ligand = bool(lig_hit and not lig_hit.get('is_metal'))
        has_metal  = bool(lig_hit and lig_hit.get('is_metal'))
        struct_score = (0.70 if has_ligand else 0.0) + (0.30 if has_metal else 0.0)
        struct_score = min(1.0, struct_score)

        # L2: database evidence
        pdbe_overlap = len(set(res_nums) & self._pdbe_res)
        pdbe_frac    = pdbe_overlap / max(len(res_nums), 1)
        db_score     = min(pdbe_frac * 4.0, 1.0)   # 25% overlap → full score

        # L3: conservation
        cs_vals   = [self._consurf[r] for r in res_nums if r in self._consurf]
        cons_mean = float(np.mean(cs_vals)) if cs_vals else -1.0
        # Scale-aware normalization across all L3 sources:
        #   HMMER (0–1 range): use directly
        #   ConSurf DB (1–9 grade): normalize to 0–1
        #   No data (-1.0): use neutral 0.5 (does not penalise)
        if cons_mean < 0.0:
            cons_score = 0.5
        elif cons_mean <= 1.0:
            cons_score = cons_mean
        else:
            cons_score = (cons_mean - 1.0) / 8.0

        nma_vals = [self._nma_fluct[r] for r in res_nums if r in self._nma_fluct]
        rcn_vals = [self._rcn_cent[r]  for r in res_nums if r in self._rcn_cent]
        nma_mean = float(np.mean(nma_vals)) if nma_vals else 0.5
        rcn_mean = float(np.mean(rcn_vals)) if rcn_vals else 0.5

        # Primary score: rigid + connected
        primary_dynamic    = 0.50 * (1.0 - nma_mean) + 0.50 * rcn_mean
        # Allosteric score: mobile + connected
        allosteric_dynamic = 0.60 * nma_mean + 0.40 * rcn_mean
        dynamic_score      = max(primary_dynamic, allosteric_dynamic)
        role_hint = ('primary'   if primary_dynamic > allosteric_dynamic
                     else 'allosteric')

        # L5: docking
        # co-crystal confirmed pockets OR rank-1 get global docking score
        has_coxtal = bool(self._closest_ligand(center, radius=6.0).get('name'))
        is_top     = pocket.get('rank', 99) == 1
        vina_score = self._best_docking_score(
            center,
            is_top_pocket=is_top,
            pocket_id=pocket.get('pocket_id'))
        dock_score = min(abs(vina_score) / 10.0, 1.0) if vina_score else 0.0

        # ── Weighted Pocket Evidence Score ────────────────────────────────
        pes = (WEIGHTS['structural']    * struct_score +
               WEIGHTS['conservation'] * cons_score   +
               WEIGHTS['dynamic']      * dynamic_score +
               WEIGHTS['database']     * db_score     +
               WEIGHTS['docking']      * dock_score)

        # Confidence label
        # Biological ligand in co-crystal → CONFIRMED regardless of PES
        lig_name   = lig_hit.get('name','') if lig_hit else ''
        is_bio_lig = lig_name.upper() in BIOLOGICAL_LIGANDS

        if has_ligand and is_bio_lig and cons_score >= 0.20:
            conf = f'🔴 CONFIRMED (co-crystal {lig_name} = known biological ligand)'
        elif has_ligand and pes >= 0.50:
            conf = f'🔴 CONFIRMED (co-crystal ligand {lig_name} detected, PES={pes:.3f})'
        elif has_metal:
            conf = '🟣 CONFIRMED (metal coordination site)'
        elif pes >= 0.55 and cons_score >= 0.60 and pdbe_overlap > 0:
            conf = '🟠 FUNCTIONALLY CONFIRMED (conservation + database)'
        # Conservation override: Cons≥0.80 + primary dynamic ≥ 0.60
        # = functionally important even without co-crystal evidence.
        # Conservation at 0.80+ means evolution has preserved these residues
        # for ~billion years — that IS biological evidence.
        elif cons_score >= 0.80 and primary_dynamic >= 0.60:
            conf = '🟠 FUNCTIONALLY SUPPORTED (highly conserved rigid hub)'
        elif pes >= 0.45 and cons_score >= 0.55:
            conf = '🟠 FUNCTIONALLY SUPPORTED (conserved lining residues)'
        elif pes >= 0.40 and role_hint == 'allosteric':
            conf = '🟡 ALLOSTERIC CANDIDATE (NMA/RCN dynamic evidence)'
        elif pes >= 0.30:
            conf = '🔵 POSSIBLE BINDING SITE (geometric + dynamic)'
        else:
            conf = '⚪ UNCERTAIN / STRUCTURAL CAVITY'

        pocket['evidence'] = {
            'PES':                  round(pes, 4),
            'confidence_label':     conf,
            'role_hint':            role_hint,
            # subscores
            'structural_score':     round(struct_score, 4),
            'conservation_score':   round(cons_score, 4),
            'consurf_mean':         round(cons_mean, 2),
            'primary_dynamic':      round(primary_dynamic, 4),
            'allosteric_dynamic':   round(allosteric_dynamic, 4),
            'dynamic_score':        round(dynamic_score, 4),
            'database_score':       round(db_score, 4),
            'docking_score':        round(dock_score, 4),
            # raw values
            'nma_mobility':         round(nma_mean, 4),
            'rcn_centrality':       round(rcn_mean, 4),
            'pdbe_overlap_n':       pdbe_overlap,
            'best_vina':            round(vina_score, 3) if vina_score else None,
            'cocrystal_hit':        lig_hit.get('name') if lig_hit else None,
            'cocrystal_dist_A':     lig_hit.get('dist') if lig_hit else None,
            'cocrystal_source':     lig_hit.get('source') if lig_hit else None,
            'is_metal_site':        has_metal,
            'lining_residues_n':    len(res_nums),
            'weights':              WEIGHTS,
        }

    def _best_docking_score(self, center: list, tol: float = 5.0,
                             is_top_pocket: bool = False,
                             pocket_id=None) -> float:
        """
        Match docking results to pocket by center distance.
        When no pocket-level coordinates available (global CSV row),
        assign the global best score to the top-ranked pocket only.
        """
        best = 0.0
        has_located = False
        for dr in self._docking_results:
            # PATH C: exact pocket_id match — highest priority
            if pocket_id is not None:
                dr_pid = dr.get('pocket_id')
                if dr_pid is not None:
                    try:
                        if int(dr_pid) == int(pocket_id):
                            for _k in ('score','best_vina','best_score','vina_score'):
                                _v = dr.get(_k)
                                if _v and float(_v) != 0.0:
                                    s = float(_v)
                                    if s < best: best = s
                                    break
                            continue
                    except (ValueError, TypeError):
                        pass
            dc = dr.get('pocket_center') or dr.get('center')
            if dc:
                d = math.sqrt(sum((a-b)**2 for a, b in zip(center, dc)))
                if d < tol:
                    s = float(dr.get('score', dr.get('best_score', 0)))
                    if s < best:
                        best = s
                    has_located = True
            elif dr.get('is_global') and is_top_pocket:
                # Global score (batch CSV, no pocket coords) → assign to top pocket
                s = float(dr.get('score', 0))
                if s < best:
                    best = s
        return best

    # ── Reclassify ──────────────────────────────────────────────────────────

    def _merge_subpockets(self, pockets: list,
                          dist_thresh: float = 10.0) -> list:
        """
        Merge pocket pairs that are physically overlapping sub-regions of
        the same binding groove. Criteria for merging:
          a) Both share the same co-crystal hit name (same ligand)
          b) Their centres are within dist_thresh Å of each other
        
        The higher-PES pocket is kept as the representative.
        The lower-PES pocket is annotated as a sub-pocket and removed
        from the ranked list, preventing false duplicate confirmations.
        
        Scientific rationale: pocket detection algorithms (fpocket,
        PolarPocket, SiteMap) often split large binding grooves into 2–3
        overlapping sub-pockets. Reporting both as independent 'CONFIRMED'
        sites overstates the evidence and misleads docking grid placement.
        """
        import math
        merged_ids = set()
        for i, pi in enumerate(pockets):
            if pi['pocket_id'] in merged_ids:
                continue
            ei  = pi.get('evidence', {})
            cxi = ei.get('cocrystal_hit')
            if not cxi:
                continue
            for j, pj in enumerate(pockets):
                if i == j or pj['pocket_id'] in merged_ids:
                    continue
                ej  = pj.get('evidence', {})
                cxj = ej.get('cocrystal_hit')
                if cxi != cxj:
                    continue
                ci, cj = pi['center'], pj['center']
                d = math.sqrt(sum((a-b)**2 for a, b in zip(ci, cj)))
                if d > dist_thresh:
                    continue
                # Merge: keep higher PES, annotate lower as sub-pocket
                keep, drop = (pi, pj) if ei.get('PES',0) >= ej.get('PES',0)                              else (pj, pi)
                drop_id = drop['pocket_id']
                keep_id = keep['pocket_id']
                merged_ids.add(drop_id)
                keep['evidence']['merged_subpocket_id'] = drop_id
                keep['evidence']['merged_dist_A']       = round(d, 2)
                keep['evidence']['confidence_label']   +=                     f' [+sub-pocket {drop_id} merged, Δ={d:.1f}Å]'
                self.log(f'Sub-pocket merge: pocket {drop_id} → '
                         f'pocket {keep_id} (Δ={d:.1f}Å, ligand={cxi})')

        # Return only non-merged pockets, re-ranked
        final = [p for p in pockets if p['pocket_id'] not in merged_ids]
        for i, p in enumerate(final):
            p['rank'] = i + 1
        return final
    def _reclassify(self, pockets: list):
        for p in pockets:
            ev   = p.get('evidence', {})
            pt   = p.get('pocket_type', {})
            code = pt.get('code', 'X')
            pes  = ev.get('PES', 0)

            if ev.get('is_metal_site'):
                pt.update({'code':'D','label':'Ion/Metal binding site',
                           'confidence':'CONFIRMED','evidence_upgraded':True})
            elif ev.get('cocrystal_hit') and not ev.get('is_metal_site'):
                pt.update({'code':'A','label':'Primary binding / active site',
                           'confidence':'CONFIRMED (co-crystal)',
                           'evidence_upgraded':True})
            elif code in ('A*','B') and pes >= 0.50 and \
                    ev.get('conservation_score',0) >= 0.60:
                pt.update({'code':'A',
                           'confidence':'FUNCTIONALLY SUPPORTED',
                           'evidence_upgraded':True})
            elif code == 'B' and ev.get('allosteric_dynamic',0) > 0.65:
                pt.update({'confidence':'ALLOSTERIC CONFIRMED (NMA+RCN)',
                           'evidence_upgraded':True})
            elif code == 'A*' and pes < 0.28:
                pt.update({'code':'X','confidence':'Low',
                           'evidence_upgraded':True})
