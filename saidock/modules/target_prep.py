import os
import subprocess
import urllib.request
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select

RCSB_URL = 'https://files.rcsb.org/download/{}.pdb'


class ProteinSelect(Select):
    def __init__(self, chain='A'):
        self.chain = chain
    def accept_model(self, model):   return 1 if model.id == 0 else 0
    def accept_chain(self, chain):   return 1 if chain.id == self.chain else 0
    def accept_residue(self, res):   return 1 if res.id[0] == ' ' else 0
    def accept_atom(self, atom):
        if atom.is_disordered() and atom.get_altloc() not in ('A', ' ', ''):
            return 0
        return 1


class TargetPrep:
    def __init__(self, target, chain='A', outdir='.', logger=None):
        self.target     = target.strip()
        self.chain      = chain.upper()
        self.outdir     = Path(outdir)
        self.logger     = logger
        self.target_id  = (target.upper()
                           if (len(target) == 4 and target.isalnum())
                           else Path(target).stem)
        self.n_residues = 0
        self._clean_pdb = None
        self._pdbqt     = None

    def log(self, msg, level='INFO'):
        if self.logger:
            self.logger.log(msg, level)

    def prepare(self) -> Path:
        self.outdir.mkdir(parents=True, exist_ok=True)
        raw   = self._download_or_copy()
        clean = self._clean(raw)
        fixed = self._fix_residues(clean)
        self._count_residues(fixed)
        self._clean_pdb = fixed
        self.log(f'Target prepared: {self.target_id}  {self.n_residues} residues', 'OK')
        return fixed

    def convert_to_pdbqt(self) -> Path:
        """
        Convert receptor PDB to PDBQT for AutoDock Vina.

        CRITICAL: use -xr flag (rigid receptor mode).
        Without -xr, obabel adds ROOT/BRANCH/ENDBRANCH tags
        which cause Vina to throw:
          "PDBQT parsing error: Unknown or inappropriate tag > ROOT"

        obabel receptor flags:
          -p 7.4     : protonate at pH 7.4
          -xr        : rigid receptor (no torsion tree tags)
          -h         : add hydrogens
          --partialcharge gasteiger : assign Gasteiger partial charges
        """
        out = self.outdir / f'{self.target_id}_receptor.pdbqt'

        # Method 1: obabel with -xr rigid flag (correct for Vina receptors)
        ret = subprocess.run([
            'obabel', str(self._clean_pdb),
            '-O', str(out),
            '-p', '7.4',
            '-xr',
            '-h',
            '--partialcharge', 'gasteiger',
            '--errorlevel', '1',
        ], capture_output=True, text=True)

        if out.exists() and out.stat().st_size > 500:
            # Verify no ROOT tags crept in
            content = out.read_text()
            if 'ROOT' not in content and 'BRANCH' not in content:
                self._pdbqt = out
                self.log(f'PDBQT generated (obabel -xr): {out.name}', 'OK')
                return out
            else:
                self.log('obabel added ROOT tags despite -xr — stripping manually', 'WARN')
                self._strip_torsion_tags(out)
                self._pdbqt = out
                self.log(f'PDBQT cleaned (stripped ROOT/BRANCH): {out.name}', 'OK')
                return out

        # Method 2: obabel without -xr, then strip tags manually
        self.log('obabel -xr failed. Trying standard conversion + manual strip.', 'WARN')
        ret2 = subprocess.run([
            'obabel', str(self._clean_pdb),
            '-O', str(out),
            '-p', '7.4',
            '-h',
            '--partialcharge', 'gasteiger',
            '--errorlevel', '1',
        ], capture_output=True, text=True)

        if out.exists() and out.stat().st_size > 500:
            self._strip_torsion_tags(out)
            self._pdbqt = out
            self.log(f'PDBQT generated + cleaned: {out.name}', 'OK')
            return out

        # Method 3: pure Python fallback — writes minimal valid receptor PDBQT
        self.log('obabel failed. Using Python PDB -> PDBQT converter.', 'WARN')
        self._python_pdb_to_pdbqt(self._clean_pdb, out)
        self._pdbqt = out
        self.log(f'PDBQT generated (Python fallback): {out.name}', 'OK')
        return out

    @staticmethod
    def _strip_torsion_tags(pdbqt_path: Path):
        """
        Remove ROOT, ENDROOT, BRANCH, ENDBRANCH lines from a receptor PDBQT.
        These tags are valid for ligands but invalid for rigid receptors.
        """
        lines = pdbqt_path.read_text().splitlines()
        clean = []
        for line in lines:
            tag = line[:6].strip()
            if tag not in ('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF'):
                clean.append(line)
        pdbqt_path.write_text('\n'.join(clean) + '\n')

    @staticmethod
    def _python_pdb_to_pdbqt(pdb_path: Path, out_path: Path):
        """
        Minimal PDB -> PDBQT conversion in pure Python.
        Copies ATOM/HETATM lines, appends atom type and charge columns.
        Sufficient for AutoDock Vina rigid receptor.
        No ROOT/BRANCH tags added.
        """
        ELEMENT_TO_VINA_TYPE = {
            'C': 'C',  'N': 'NA', 'O': 'OA', 'S': 'SA',
            'P': 'P',  'H': 'HD', 'F': 'F',  'CL': 'CL',
            'BR': 'BR','I': 'I',  'FE': 'FE','ZN': 'ZN',
            'CA': 'CA','MG': 'MG',
        }
        out_lines = []
        with open(pdb_path) as f:
            for line in f:
                record = line[:6].strip()
                if record not in ('ATOM', 'HETATM'):
                    continue
                # Pad line to 80 chars
                line = line.rstrip('\n').ljust(80)
                elem = line[76:78].strip().upper()
                if not elem:
                    name = line[12:16].strip().lstrip('0123456789')
                    elem = name[:2].upper() if len(name) >= 2 else name.upper()
                vtype = ELEMENT_TO_VINA_TYPE.get(elem, 'A')
                charge = '0.000'
                # Write PDBQT line: original 66 chars + charge + type
                pdbqt_line = f'{line[:66]}  {float(charge):6.3f} {vtype:<2}'
                out_lines.append(pdbqt_line)
        out_path.write_text('\n'.join(out_lines) + '\n')

    def _download_or_copy(self) -> Path:
        out = self.outdir / f'{self.target_id}_raw.pdb'
        if Path(self.target).is_file():
            import shutil
            shutil.copy(self.target, out)
            self.log(f'Copied local PDB: {self.target}')
        else:
            if not out.exists():
                url = RCSB_URL.format(self.target_id)
                self.log(f'Downloading from RCSB: {url}')
                try:
                    urllib.request.urlretrieve(url, out)
                except Exception as e:
                    raise RuntimeError(f'PDB download failed ({url}): {e}')
            else:
                self.log(f'Using cached PDB: {out.name}')
        return out

    def _clean(self, pdb_path: Path) -> Path:
        out = self.outdir / f'{self.target_id}_chain{self.chain}.pdb'
        try:
            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(self.target_id, str(pdb_path))
            io = PDBIO()
            io.set_structure(struct)
            io.save(str(out), ProteinSelect(self.chain))
            self.log(f'Chain {self.chain} extracted: {out.name}')
        except Exception as e:
            self.log(f'Chain extraction failed: {e}', 'WARN')
            import shutil
            shutil.copy(pdb_path, out)
        return out

    def _fix_residues(self, pdb_path: Path) -> Path:
        out = self.outdir / f'{self.target_id}_fixed.pdb'
        try:
            import pdbfixer
            from openmm.app import PDBFile
            fixer = pdbfixer.PDBFixer(filename=str(pdb_path))
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.removeHeterogens(True)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            with open(out, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
            self.log(f'Missing residues repaired: {out.name}')
        except ImportError:
            self.log('pdbfixer not installed — skipping repair', 'WARN')
            import shutil
            shutil.copy(pdb_path, out)
        except Exception as e:
            self.log(f'pdbfixer failed: {e}', 'WARN')
            import shutil
            shutil.copy(pdb_path, out)
        return out

    def _count_residues(self, pdb_path: Path):
        try:
            parser = PDBParser(QUIET=True)
            struct = parser.get_structure('X', str(pdb_path))
            self.n_residues = len(list(struct.get_residues()))
        except Exception:
            self.n_residues = 0
