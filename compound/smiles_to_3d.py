from rdkit import Chem
from rdkit.Chem import AllChem

# Flavopiridol SMILES
smiles = "CC(C)NCCNC(=O)C1=CC=C(C=C1)O"

# Convert SMILES into molecule
mol = Chem.MolFromSmiles(smiles)

if mol is None:
    print("Invalid SMILES")
    exit()

# Add hydrogen atoms
mol = Chem.AddHs(mol)

# Generate 3D coordinates
result = AllChem.EmbedMolecule(mol, randomSeed=42)

if result != 0:
    print("Could not generate 3D structure")
    exit()

# Optimize the 3D structure
AllChem.UFFOptimizeMolecule(mol)

# Save the molecule as SDF
writer = Chem.SDWriter("compound/flavopiridol_3d.sdf")
writer.write(mol)
writer.close()

print("3D structure generated successfully.")
print("Saved as: compound/flavopiridol_3d.sdf")