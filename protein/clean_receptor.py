from Bio.PDB import PDBParser, PDBIO, Select

input_pdb = "protein/CDK1_chainA.pdb"
output_pdb = "protein/CDK1_clean.pdb"


class ProteinOnlySelect(Select):

    def accept_residue(self, residue):
        # Keep only standard amino-acid residues
        hetflag = residue.id[0]

        if hetflag == " ":
            return 1

        return 0


parser = PDBParser(QUIET=True)
structure = parser.get_structure("CDK1", input_pdb)

io = PDBIO()
io.set_structure(structure)
io.save(output_pdb, ProteinOnlySelect())

print("Clean protein receptor created successfully.")
print("Saved as:", output_pdb)