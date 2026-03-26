"""
Stage 05: Protein-Ligand Interaction Fingerprinting
Primary:  ProLIF 2.x (MDAnalysis-based, corrected API for v2.1+)
Fallback: Distance-based H-bond + hydrophobic contact detection
"""
import warnings
warnings.filterwarnings(
    'ignore',
    message='MDAnalysis.topology.tables has been moved',
    category=DeprecationWarning,
)
warnings.filterwarnings(
    'ignore',
    message='please use MorganGenerator',
    category=DeprecationWarning,
)
import json
import numpy as np
from pathlib import Path



def _prolif_worker(receptor_pdb, ligand_pdbqt, outdir):
    """
    Standalone ProLIF worker — runs in separate process for timeout safety.
    Returns dict with hbonds, hydrophobic, interactions, method.
    """
    import subprocess, tempfile
    from pathlib import Path
    import MDAnalysis as mda
    import prolif as plf

    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        lig_pdb = tmp.name
    with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False) as tmp2:
        lig_sdf = tmp2.name

    subprocess.run(['obabel', ligand_pdbqt, '-O', lig_sdf, '--errorlevel', '1'],
                   capture_output=True, text=True)
    # Use RDKit for SDF→PDB (avoids valence errors from obabel -h on PDBQT-derived SDF)
    try:
        from rdkit import Chem
        _mol = Chem.MolFromMolFile(lig_sdf, sanitize=False, removeHs=False)
        if _mol is None:
            raise ValueError('RDKit could not read ligand SDF')
        try:
            Chem.SanitizeMol(_mol)
        except Exception:
            Chem.SanitizeMol(
                _mol,
                Chem.SanitizeFlags.SANITIZE_ALL ^
                Chem.SanitizeFlags.SANITIZE_PROPERTIES
            )
        Chem.MolToPDBFile(_mol, lig_pdb)
    except Exception as _rdkit_err:
        # Final fallback: obabel without -h flag
        subprocess.run(
            ['obabel', lig_sdf, '-O', lig_pdb, '--errorlevel', '1'],
            capture_output=True, text=True
        )
    if not Path(lig_pdb).exists():
        raise FileNotFoundError(f"Ligand PDB not created: {lig_pdb}")

    rec_lines, lig_lines = [], []
    with open(receptor_pdb) as f:
        for line in f:
            if line[:4] in ('ATOM', 'HETA', 'TER ', 'END '):
                rec_lines.append(line)
    with open(lig_pdb) as f:
        for line in f:
            if line[:4] in ('ATOM', 'HETA'):
                lig_lines.append(line[:17] + 'LIG' + line[20:])

    complex_pdb = str(Path(outdir) / 'complex_for_prolif.pdb')
    with open(complex_pdb, 'w') as f:
        f.writelines(rec_lines)
        f.write('TER\n')
        f.writelines(lig_lines)
        f.write('END\n')

    u    = mda.Universe(complex_pdb)
    lig  = u.select_atoms('resname LIG')
    prot = u.select_atoms('protein')

    if len(lig) == 0 or len(prot) == 0:
        raise ValueError("Empty ligand or protein selection")

    fp = plf.Fingerprint(['HBDonor', 'HBAcceptor', 'Hydrophobic',
                          'PiStacking', 'Cationic', 'Anionic'])
    try:
        fp.run(u.trajectory, lig, prot, n_jobs=1)
    except TypeError:
        fp.run(lig, prot)

    df = fp.to_dataframe()
    if df.empty:
        return {'hbonds': 0, 'hydrophobic': 0, 'interactions': [], 'method': 'prolif_2x_no_contacts'}

    hbonds, hydrophobic, contacts = 0, 0, []
    for col in df.columns:
        if df[col].any():
            residue    = str(col[0]) if isinstance(col, tuple) else str(col)
            inter_type = str(col[1]) if isinstance(col, tuple) and len(col) > 1 else 'unknown'
            contacts.append({'residue': residue, 'type': inter_type})
            if 'HB' in inter_type.upper():
                hbonds += 1
            if 'HYDRO' in inter_type.upper() or 'PI' in inter_type.upper():
                hydrophobic += 1

    return {'hbonds': hbonds, 'hydrophobic': hydrophobic,
            'interactions': contacts, 'method': 'prolif_2x'}


class InteractionAnalyzer:
    def __init__(self, receptor_pdb, ligand_pose_pdbqt=None,
                 outdir='.', logger=None, docking_results=None,
                 **kwargs):
        """
        receptor_pdb        : path to fixed receptor PDB
        ligand_pose_pdbqt   : path to best docked pose PDBQT
        docking_results     : list of docking result dicts (backward compat)
                              — best pose file extracted automatically
        """
        self.receptor_pdb = Path(receptor_pdb)
        self.outdir       = Path(outdir)
        self.logger       = logger

        # Resolve ligand pose — direct path takes priority
        if ligand_pose_pdbqt:
            self.ligand_pose_pdbqt = Path(ligand_pose_pdbqt)
        elif docking_results:
            # Extract best pose file from docking_results (old calling convention)
            valid = [r for r in docking_results
                     if r.get('best_score', 0) < 0 and r.get('pose_file')]
            if valid:
                best = min(valid, key=lambda r: r['best_score'])
                self.ligand_pose_pdbqt = Path(best['pose_file'])
            else:
                self.ligand_pose_pdbqt = None
        else:
            self.ligand_pose_pdbqt = None

    def log(self, msg, level='INFO'):
        if self.logger:
            self.logger.log(msg, level)

    def analyse(self) -> dict:
        self.outdir.mkdir(parents=True, exist_ok=True)

        if self.ligand_pose_pdbqt and self.ligand_pose_pdbqt.exists():
            result = self._try_prolif()
            if result is None:
                result = self._distance_fallback()
        else:
            self.log('No pose file — using receptor-only distance estimation', 'WARN')
            result = {'hbonds': 0, 'hydrophobic': 0,
                      'interactions': [], 'method': 'no_pose'}

        out = self.outdir / 'interactions.json'
        with open(out, 'w') as f:
            json.dump(result, f, indent=2)
        self.log(
            f'Interactions: H-bonds={result["hbonds"]}  '
            f'Hydrophobic={result["hydrophobic"]}',
            'OK'
        )
        return result

    # ------------------------------------------------------------------
    def _try_prolif(self) -> dict:
        """
        ProLIF 2.x API with 30s hard timeout via multiprocessing.
        Falls back to distance-based method on timeout or valence error.
        """
        import multiprocessing as mp

        def _run_prolif(result_queue, receptor_pdb, ligand_pdbqt, outdir):
            try:
                res = _prolif_worker(receptor_pdb, ligand_pdbqt, outdir)
                result_queue.put(res)
            except Exception as e:
                result_queue.put(e)

        q   = mp.Queue()
        proc = mp.Process(
            target=_run_prolif,
            args=(q, str(self.receptor_pdb),
                  str(self.ligand_pose_pdbqt), str(self.outdir)),
            daemon=True,
        )
        proc.start()
        proc.join(timeout=120)          # 30 seconds hard limit

        if proc.is_alive():
            proc.terminate()
            proc.join()
            self.log('ProLIF timed out after 30s — using distance-based fallback',
                     'WARN')
            return None

        if not q.empty():
            result = q.get()
            if isinstance(result, Exception):
                self.log(f'ProLIF failed: {result} — using distance-based fallback',
                         'WARN')
                return None
            return result
        return None

    def _try_prolif_inner(self) -> dict:
        """Original ProLIF logic — called inside subprocess."""
        try:
            import MDAnalysis as mda
            import prolif as plf

            # Convert PDBQT ligand pose to PDB for MDAnalysis
            import subprocess, tempfile
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
                lig_pdb = tmp.name
            # Convert PDBQT → SDF first (preserves aromaticity),
            # then SDF → PDB (avoids explicit-valence errors from direct PDBQT→PDB)
            import tempfile
            with tempfile.NamedTemporaryFile(suffix='.sdf', delete=False) as tmp2:
                lig_sdf = tmp2.name
            subprocess.run(
                ['obabel', str(self.ligand_pose_pdbqt),
                 '-O', lig_sdf, '--errorlevel', '1'],
                capture_output=True, text=True
            )
            if Path(lig_sdf).exists() and Path(lig_sdf).stat().st_size > 0:
                subprocess.run(
                    ['obabel', lig_sdf, '-O', lig_pdb,
                     '--errorlevel', '1', '-h'],
                    capture_output=True, text=True
                )
            else:
                # Fallback: direct conversion
                subprocess.run(
                    ['obabel', str(self.ligand_pose_pdbqt),
                     '-O', lig_pdb, '--errorlevel', '1'],
                    capture_output=True, text=True
                )
            if not Path(lig_pdb).exists():
                return None

            # Build complex PDB
            complex_pdb = str(self.outdir / 'complex_for_prolif.pdb')
            rec_lines   = []
            with open(self.receptor_pdb) as f:
                for line in f:
                    if line[:4] in ('ATOM', 'HETA', 'TER ', 'END '):
                        rec_lines.append(line)
            lig_lines = []
            lig_resnum = 1
            with open(lig_pdb) as f:
                for line in f:
                    if line[:4] in ('ATOM', 'HETA'):
                        # Rename residue to LIG for ProLIF
                        line = line[:17] + 'LIG' + line[20:]
                        lig_lines.append(line)
            with open(complex_pdb, 'w') as f:
                f.writelines(rec_lines)
                f.write('TER\n')
                f.writelines(lig_lines)
                f.write('END\n')

            # Sanitize ligand valence issues before handing to ProLIF
            try:
                from rdkit import Chem
                tmp_mol = Chem.MolFromPDBFile(lig_pdb, sanitize=False, removeHs=False)
                if tmp_mol is not None:
                    try:
                        Chem.SanitizeMol(tmp_mol)
                    except Exception:
                        # Remove Hs causing valence errors and re-write
                        tmp_mol = Chem.RemoveHs(tmp_mol, sanitize=False)
                        Chem.SanitizeMol(tmp_mol,
                            Chem.SanitizeFlags.SANITIZE_ALL ^
                            Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                    writer = Chem.PDBWriter(lig_pdb)
                    writer.write(tmp_mol)
                    writer.close()
            except Exception:
                pass  # proceed with original file

            # ProLIF 2.x: load as MDAnalysis universe
            u    = mda.Universe(complex_pdb)
            lig  = u.select_atoms('resname LIG')
            prot = u.select_atoms('protein')

            if len(lig) == 0 or len(prot) == 0:
                return None

            fp = plf.Fingerprint([
                'HBDonor', 'HBAcceptor', 'Hydrophobic',
                'PiStacking', 'Cationic', 'Anionic',
            ])

            # ProLIF 2.1+ API: run() takes (universe_or_trajectory, lig_ag, prot_ag)
            try:
                fp.run(u.trajectory, lig, prot, n_jobs=1)
            except TypeError:
                # Older 2.x fallback
                fp.run(lig, prot)

            df = fp.to_dataframe()
            if df.empty:
                return {'hbonds': 0, 'hydrophobic': 0,
                        'interactions': [], 'method': 'prolif_2x_no_contacts'}

            hbonds     = 0
            hydrophobic = 0
            contacts   = []
            for col in df.columns:
                if df[col].any():
                    residue    = str(col[0]) if isinstance(col, tuple) else str(col)
                    inter_type = str(col[1]) if isinstance(col, tuple) and len(col)>1 else 'unknown'
                    contacts.append({'residue': residue, 'type': inter_type})
                    if 'HB' in inter_type.upper():
                        hbonds += 1
                    if 'HYDRO' in inter_type.upper() or 'PI' in inter_type.upper():
                        hydrophobic += 1

            return {
                'hbonds':      hbonds,
                'hydrophobic': hydrophobic,
                'interactions': contacts,
                'method':      'prolif_2x',
            }

        except Exception as e:
            self.log(f'ProLIF failed: {e} — using distance-based fallback', 'WARN')
            return None

    # ------------------------------------------------------------------
    def _distance_fallback(self) -> dict:
        """
        Pure-Python distance-based interaction detection.
        Uses receptor PDB + ligand PDBQT coordinates directly.
        """
        rec_coords, rec_elems, rec_res = self._read_pdb_atoms(self.receptor_pdb)
        lig_coords, lig_elems          = self._read_pdbqt_atoms(self.ligand_pose_pdbqt)

        if len(rec_coords) == 0 or len(lig_coords) == 0:
            return {'hbonds': 0, 'hydrophobic': 0,
                    'interactions': [], 'method': 'distance_fallback_no_coords'}

        POLAR    = {'N', 'O', 'F', 'S'}
        HYDRO    = {'C', 'S'}
        HBOND_D  = 3.5   # A
        HYDRO_D  = 4.5   # A

        hb_contacts   = []
        hydr_contacts = []

        for i, (lc, le) in enumerate(zip(lig_coords, lig_elems)):
            diffs   = rec_coords - lc
            dists   = np.sqrt((diffs**2).sum(axis=1))
            close   = np.where(dists < max(HBOND_D, HYDRO_D))[0]
            for j in close:
                re = rec_elems[j].upper()
                d  = float(dists[j])
                res_label = rec_res[j]
                if le.upper() in POLAR and re in POLAR and d < HBOND_D:
                    hb_contacts.append({
                        'residue': res_label, 'type': 'HBond',
                        'distance': round(d, 2),
                    })
                if le.upper() in HYDRO and re in HYDRO and d < HYDRO_D:
                    hydr_contacts.append({
                        'residue': res_label, 'type': 'Hydrophobic',
                        'distance': round(d, 2),
                    })

        # Deduplicate by residue
        hb_unique   = list({c['residue']: c for c in hb_contacts}.values())
        hy_unique   = list({c['residue']: c for c in hydr_contacts}.values())

        return {
            'hbonds':       len(hb_unique),
            'hydrophobic':  len(hy_unique),
            'interactions': hb_unique + hy_unique,
            'method':       'distance_fallback',
        }

    @staticmethod
    def _read_pdb_atoms(pdb_path):
        coords, elems, residues = [], [], []
        with open(pdb_path) as f:
            for line in f:
                if line[:4] not in ('ATOM', 'HETA'):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    elem = line[76:78].strip().upper()
                    if not elem:
                        elem = line[12:16].strip().lstrip('0123456789')[:1]
                    res  = f"{line[17:20].strip()}{line[22:26].strip()}"
                    coords.append([x, y, z])
                    elems.append(elem)
                    residues.append(res)
                except (ValueError, IndexError):
                    continue
        return np.array(coords), elems, residues

    # AutoDock atom type → element (for H-bond polar detection)
    AUTODOCK_ELEM = {
        'HD': 'H', 'HS': 'H',               # polar H donors
        'OA': 'O', 'OS': 'O',               # H-bond acceptor O
        'NA': 'N', 'NS': 'N',               # H-bond acceptor N
        'SA': 'S',                           # acceptor S
        'C' : 'C', 'A' : 'C',               # aliphatic / aromatic C
        'N' : 'N', 'O' : 'O', 'S' : 'S',
        'H' : 'H', 'F' : 'F', 'CL': 'CL',
        'BR': 'BR', 'I' : 'I', 'FE': 'FE',
        'ZN': 'ZN', 'MG': 'MG', 'MN': 'MN',
        'P' : 'P',
    }

    @classmethod
    def _read_pdbqt_atoms(cls, pdbqt_path):
        coords, elems = [], []
        with open(pdbqt_path) as f:
            for line in f:
                if line[:4] not in ('ATOM', 'HETA'):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    # AutoDock atom type is in cols 77-79 (1-indexed)
                    ad_type = line[77:79].strip().upper()
                    elem    = cls.AUTODOCK_ELEM.get(ad_type, ad_type[:1])
                    if not elem:
                        elem = line[12:16].strip().lstrip('0123456789')[:1].upper()
                    coords.append([x, y, z])
                    elems.append(elem)
                except (ValueError, IndexError):
                    continue
        return np.array(coords), elems
