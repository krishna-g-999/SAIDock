import urllib.request
import urllib.parse
import json

PUBCHEM_URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'

def name_to_smiles(name: str) -> dict:
    """
    Fetch canonical SMILES, InChIKey, MW, and CID from PubChem by compound name.
    Returns dict with keys: smiles, inchikey, mw, cid, iupac_name
    """
    encoded = urllib.parse.quote(name)
    url     = f'{PUBCHEM_URL}/name/{encoded}/property/CanonicalSMILES,InChIKey,MolecularWeight,IUPACName/JSON'
    try:
        with urllib.request.urlopen(url, timeout=15) as r:
            data = json.loads(r.read())
        props = data['PropertyTable']['Properties'][0]
        return {
            'smiles':      props.get('CanonicalSMILES', ''),
            'inchikey':    props.get('InChIKey', ''),
            'mw':          props.get('MolecularWeight', 0),
            'cid':         props.get('CID', ''),
            'iupac_name':  props.get('IUPACName', name),
        }
    except Exception as e:
        raise RuntimeError(f'PubChem lookup failed for "{name}": {e}')

def smiles_to_2d_image_url(smiles: str, size: int = 300) -> str:
    encoded = urllib.parse.quote(smiles)
    return f'{PUBCHEM_URL}/smiles/{encoded}/PNG?image_size={size}x{size}'
