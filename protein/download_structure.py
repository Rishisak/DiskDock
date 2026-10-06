import requests

pdb_id = "6GU2"

url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

response = requests.get(url)

if response.status_code == 200:
    with open(f"protein/{pdb_id}.pdb", "w") as file:
        file.write(response.text)

    print(f"{pdb_id}.pdb downloaded successfully.")

else:
    print("Failed to download PDB structure.")