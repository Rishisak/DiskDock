import requests

def get_smiles_from_pubchem(compound_name):

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        f"compound/name/{compound_name}/property/"
        "CanonicalSMILES/JSON"
    )

    response = requests.get(url)

    if response.status_code != 200:
        raise ValueError(
            f"Could not find compound: {compound_name}"
        )

    data = response.json()

    smiles = data["PropertyTable"]["Properties"][0]["ConnectivitySMILES"]

    return smiles