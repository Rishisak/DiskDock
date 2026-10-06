import requests

uniprot_id = "P06493"

url = (
    "https://search.rcsb.org/rcsbsearch/v2/query"
)

query = {
    "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
            "operator": "exact_match",
            "value": "UniProt"
        }
    },
    "return_type": "entry",
    "request_options": {
        "paginate": {
            "start": 0,
            "rows": 100
        }
    }
}

response = requests.post(url, json=query)

if response.status_code != 200:
    print("Error searching RCSB PDB")
    print(response.text)
    exit()

data = response.json()

print("PDB structures found:")

for result in data.get("result_set", []):
    print(result["identifier"])