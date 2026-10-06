from compound.structure import smiles_to_3d


smiles = "CN1CCC(C(C1)O)C2=C(C=C(C3=C2OC(=CC3=O)C4=CC=CC=C4Cl)O)O"

smiles_to_3d(
    smiles,
    "results/test_flavopiridol.sdf"
)