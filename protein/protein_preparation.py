from Bio.PDB import PDBParser, PDBIO, Select


class ChainSelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"


# Input PDB
input_pdb = "protein/6GU2.pdb"

# Output PDB containing only CDK1 chain A
output_pdb = "protein/CDK1_chainA.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("CDK1", input_pdb)

io = PDBIO()
io.set_structure(structure)
io.save(output_pdb, ChainSelect())

print("CDK1 Chain A extracted successfully.")
print("Saved as:", output_pdb)