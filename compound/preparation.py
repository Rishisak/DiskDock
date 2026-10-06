from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy


def prepare_ligand(input_sdf, output_pdbqt):
    """
    Convert a 3D SDF ligand into PDBQT format using Meeko.
    """

    print("\n[5] Loading ligand...")

    supplier = Chem.SDMolSupplier(
        input_sdf,
        removeHs=False
    )

    mol = supplier[0]

    if mol is None:
        raise ValueError("Could not load ligand from SDF.")

    print("Ligand loaded successfully.")

    print("[6] Preparing ligand with Meeko...")

    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)

    pdbqt_string, is_ok, error_msg = (
        PDBQTWriterLegacy.write_string(mol_setups[0])
    )

    if not is_ok:
        raise RuntimeError(
            f"Meeko ligand preparation failed: {error_msg}"
        )

    with open(output_pdbqt, "w") as f:
        f.write(pdbqt_string)

    print("Ligand preparation completed.")
    print("Saved as:", output_pdbqt)

    return output_pdbqt