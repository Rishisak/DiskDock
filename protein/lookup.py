import requests
import os


def search_pdb(protein_name):
    """
    Search RCSB PDB for structures matching the protein name.
    """

    print("\n[Protein Lookup]")
    print("Searching RCSB PDB for:", protein_name)

    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": protein_name
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 10
            }
        }
    }

    response = requests.post(
        url,
        json=query,
        timeout=30
    )

    if response.status_code != 200:
        raise RuntimeError(
            f"RCSB PDB search failed. "
            f"Status code: {response.status_code}"
        )

    data = response.json()

    results = data.get("result_set", [])

    if not results:
        raise ValueError(
            f"No PDB structures found for protein: {protein_name}"
        )

    print("\nStructures found:")

    for i, result in enumerate(results, start=1):

        pdb_id = result["identifier"]
        score = result.get("score", 0)

        print(
            f"{i}. {pdb_id}  score={score}"
        )

    # -------------------------------------------------
    # Structure selection
    # -------------------------------------------------

    pdb_ids = [
        result["identifier"]
        for result in results
    ]

    # Temporary preferred structure for CDK1.
    # 6GU2 contains the co-crystallized F9Z ligand
    # used in our CDK1-Flavopiridol test.
    preferred_pdb = "6GU2"

    if (
        protein_name.lower() == "cdk1"
        and preferred_pdb in pdb_ids
    ):
        pdb_id = preferred_pdb

        print(
            "\nPreferred CDK1 structure found."
        )

    else:
        pdb_id = results[0]["identifier"]

    print("\nSelected PDB:", pdb_id)

    return pdb_id


def download_pdb(pdb_id, output_file):
    """
    Download a PDB structure from RCSB.
    """

    print("\nDownloading PDB structure...")

    url = (
        f"https://files.rcsb.org/download/"
        f"{pdb_id}.pdb"
    )

    response = requests.get(
        url,
        timeout=30
    )

    if response.status_code != 200:
        raise RuntimeError(
            f"Could not download PDB structure: {pdb_id}"
        )

    # Create folder if necessary
    directory = os.path.dirname(output_file)

    if directory:
        os.makedirs(
            directory,
            exist_ok=True
        )

    with open(
        output_file,
        "w"
    ) as file:

        file.write(response.text)

    print(
        "PDB downloaded successfully."
    )

    print(
        "Saved as:",
        output_file
    )

    return output_file


def get_protein_structure(protein_name):
    """
    Search for a protein structure and download it.
    """

    pdb_id = search_pdb(
        protein_name
    )

    output_file = (
        f"protein/{pdb_id}.pdb"
    )

    download_pdb(
        pdb_id,
        output_file
    )

    return output_file, pdb_id