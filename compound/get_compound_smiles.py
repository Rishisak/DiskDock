import requests

compound_name = "Flavopiridol"

url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/CanonicalSMILES/JSON"

response = requests.get(url)

if response.status_code == 200:
    data = response.json()

    smiles = data["PropertyTable"]["Properties"][0]["ConnectivitySMILES"]

    print("Compound:", compound_name)
    print("SMILES:", smiles)

else:
    print("Could not find the compound.")