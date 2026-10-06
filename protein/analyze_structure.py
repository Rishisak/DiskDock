from Bio.PDB import PDBParser

pdb_file = "protein/6GU2.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("CDK1", pdb_file)

print("===== 6GU2 STRUCTURE =====")

for model in structure:
    print(f"\nModel: {model.id}")

    for chain in model:
        protein_residues = []

        for residue in chain:
            # Ignore water and other HETATM residues
            if residue.id[0] == " ":
                protein_residues.append(residue)

        if protein_residues:
            first = protein_residues[0]
            last = protein_residues[-1]

            print(
                f"Chain {chain.id}: "
                f"{len(protein_residues)} residues | "
                f"{first.resname}{first.id[1]} -> "
                f"{last.resname}{last.id[1]}"
            )

print("\n===== END =====")