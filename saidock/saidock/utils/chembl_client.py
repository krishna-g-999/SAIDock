import urllib.request
import urllib.parse
import json

CHEMBL_BASE  = 'https://www.ebi.ac.uk/chembl/api/data'
UNIPROT_BASE = 'https://www.ebi.ac.uk/pdbe/api/mappings/uniprot'

# Known PDB -> UniProt -> ChEMBL mappings (fast lookup without API round-trip)
PDB_CHEMBL_MAP = {
    # ── CK2 alpha (CSNK2A1) ─────────────────────────────────────────
    '3PE1': 'CHEMBL2185', '2PVR': 'CHEMBL2185', '3NSZ': 'CHEMBL2185',
    '3KXM': 'CHEMBL2185', '1J3H': 'CHEMBL2185',
    # ── CDK6 ────────────────────────────────────────────────────────
    '5L2I': 'CHEMBL2828', '4L7B': 'CHEMBL2828', '2HBJ': 'CHEMBL2828',
    '2EUF': 'CHEMBL2828', '1XO2': 'CHEMBL2828',
    # ── KEAP1 / NRF2 ────────────────────────────────────────────────
    '4XMB': 'CHEMBL4523', '3EQC': 'CHEMBL4523', '5FNQ': 'CHEMBL4523',
    '6QIB': 'CHEMBL4523', '5CGJ': 'CHEMBL4523', '4ZY3': 'CHEMBL4523',
    # ── TDP-43 RRM1 (TARDBP) — ChEMBL4523419 ───────────────────────
    # Note: TDP-43 has <30 direct small-molecule entries in ChEMBL.
    # Model is trained on RNA-recognition domain binders + decoys.
    '5HDN': 'CHEMBL4523419', '2N3X': 'CHEMBL4523419',
    '4IUF': 'CHEMBL4523419', '4BS2': 'CHEMBL4523419',
    '7XS0': 'CHEMBL4523419', '2KFN': 'CHEMBL4523419',
    # ── TDP-43 C-terminal / aggregation peptide ──────────────────────
    '6N37': 'CHEMBL4523419', '5WHN': 'CHEMBL4523419',
    '6CFH': 'CHEMBL4523419',
    # ── Parkin / PARK2 (UniProt O60260) ─────────────────────────────
    '2C9V': 'CHEMBL3471',   '5C1Z': 'CHEMBL3471',
    '4K7D': 'CHEMBL3471',   '4BM9': 'CHEMBL3471',
    # ── Huntingtin / HTT ────────────────────────────────────────────
    '6EZ8': 'CHEMBL4630',   '3IO4': 'CHEMBL4630',
    # ── SOD1 (ALS) ──────────────────────────────────────────────────
    '1PU0': 'CHEMBL2581',   '2C9V': 'CHEMBL2581',
    # ── Negative control ─────────────────────────────────────────────
    '1AKI': '',
}


def pdb_to_chembl_id(pdb_id: str) -> str:
    """
    Convert PDB ID to ChEMBL target ID.
    Strategy: check local map first, then try PDB -> UniProt -> ChEMBL API chain.
    """
    pdb_id = pdb_id.upper().strip()

    # Step 1: local map (fast, no API calls)
    if pdb_id in PDB_CHEMBL_MAP and PDB_CHEMBL_MAP[pdb_id]:
        return PDB_CHEMBL_MAP[pdb_id]

    # Step 2: PDB -> UniProt via PDBe API
    try:
        url = f'{UNIPROT_BASE}/{pdb_id.lower()}'
        with urllib.request.urlopen(url, timeout=10) as r:
            data   = json.loads(r.read())
        uniprot_ids = list(data.get(pdb_id.lower(), {}).keys())
        if not uniprot_ids:
            return ''
        uniprot_id  = uniprot_ids[0]
    except Exception:
        return ''

    # Step 3: UniProt -> ChEMBL
    return get_target_chembl_id(uniprot_id)


def get_target_chembl_id(uniprot_id: str) -> str:
    url = (f'{CHEMBL_BASE}/target?'
           f'target_components__accession={uniprot_id}&format=json')
    try:
        with urllib.request.urlopen(url, timeout=10) as r:
            data = json.loads(r.read())
        targets = data.get('targets', [])
        for t in targets:
            if t.get('target_type') == 'SINGLE PROTEIN':
                return t.get('target_chembl_id', '')
        if targets:
            return targets[0].get('target_chembl_id', '')
    except Exception:
        pass
    return ''


def fetch_target_activities(target_id: str,
                             assay_types: list = None,
                             max_records: int = 200) -> list:
    """
    Fetch experimental binding activities from ChEMBL.
    Accepts PDB ID or ChEMBL ID directly.
    """
    if assay_types is None:
        assay_types = ['IC50', 'Kd', 'Ki', 'EC50']

    # Resolve to ChEMBL ID if a PDB ID was given
    chembl_id = target_id
    if not target_id.startswith('CHEMBL'):
        chembl_id = pdb_to_chembl_id(target_id)
    if not chembl_id:
        return []

    types_str = ','.join(assay_types)
    url = (f'{CHEMBL_BASE}/activity?'
           f'target_chembl_id={chembl_id}'
           f'&standard_type__in={types_str}'
           f'&standard_units=nM'
           f'&pchembl_value__isnull=false'
           f'&limit={max_records}&format=json')
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            data = json.loads(r.read())
        results = []
        for act in data.get('activities', []):
            smiles  = act.get('canonical_smiles', '')
            pchembl = act.get('pchembl_value')
            stdval  = act.get('standard_value')
            if smiles and pchembl:
                try:
                    results.append({
                        'smiles':         smiles,
                        'pchembl_value':  float(pchembl),
                        'standard_value': float(stdval) if stdval else None,
                        'dg_estimated':   round(-0.592 * float(pchembl), 3),
                        'assay_type':     act.get('standard_type', ''),
                        'chembl_id':      chembl_id,
                    })
                except (ValueError, TypeError):
                    continue
        return results
    except Exception:
        return []
