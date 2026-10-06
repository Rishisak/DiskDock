from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

input_file = "compound/flavopiridol_3d.sdf"
output_file = "compound/flavopiridol.pdbqt"

# Read ligand
supplier = Chem.SDMolSupplier(input_file, removeHs=False)
mol = supplier[0]

if mol is None:
    print("Could not read ligand.")
    exit()

print("Ligand loaded successfully.")

# Prepare ligand
preparator = MoleculePreparation()
mol_setups = preparator.prepare(mol)

# Convert prepared molecule to PDBQT
pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(
    mol_setups[0]
)

if not is_ok:
    print("PDBQT conversion failed:")
    print(error_msg)
    exit()

# Save PDBQT
with open(output_file, "w") as f:
    f.write(pdbqt_string)

print("Ligand preparation completed.")
print("Saved as:", output_file)