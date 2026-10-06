from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_3d(smiles, output_file):
    """
    Convert a SMILES string into a 3D SDF structure.

    Parameters:
        smiles       : input compound SMILES
        output_file  : path where the SDF file will be saved
    """

    print("\n[1] Converting SMILES to molecule...")

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Invalid SMILES string.")

    print("Molecule created successfully.")

    # Add hydrogens
    mol = Chem.AddHs(mol)

    print("[2] Generating 3D coordinates...")

    # Generate 3D coordinates
    result = AllChem.EmbedMolecule(
        mol,
        randomSeed=42
    )

    if result != 0:
        raise RuntimeError("Could not generate 3D coordinates.")

    print("3D coordinates generated.")

    # Optimize geometry
    print("[3] Optimizing molecular geometry...")

    result = AllChem.UFFOptimizeMolecule(mol)

    if result != 0:
        print("Warning: UFF optimization did not fully converge.")

    print("Geometry optimization completed.")

    # Save SDF
    writer = Chem.SDWriter(str(output_file))
    writer.write(mol)
    writer.close()

    print("[4] 3D structure saved:")
    print(output_file)

    return output_file
