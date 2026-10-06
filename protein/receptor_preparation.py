from meeko import MoleculePreparation
from rdkit import Chem


input_file = "protein/CDK1_chainA.pdb"
output_file = "protein/CDK1_receptor.pdbqt"


# Read protein
mol = Chem.MolFromPDBFile(
    input_file,
    removeHs=False,
    sanitize=True
)

if mol is None:
    print("Could not read protein structure.")
    exit()

print("Protein loaded successfully.")


# Add explicit hydrogens
mol = Chem.AddHs(mol)

print("Explicit hydrogens added.")


# Prepare receptor
preparator = MoleculePreparation()

try:
    setups = preparator.prepare(mol)
except Exception as e:
    print("Protein preparation failed:")
    print(e)
    exit()


# Write PDBQT
setup = setups[0]

pdbqt_string = setup.write_pdbqt_string()


with open(output_file, "w") as f:
    f.write(pdbqt_string)


print("Protein preparation completed.")
print("Saved as:", output_file)