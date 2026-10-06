from docking.box import calculate_box_from_ligand

box = calculate_box_from_ligand(
    "protein/6GU2.pdb",
    "F9Z",
    chain_id="A",
    padding=10.0
)

print("\nFinal box:")
print(box)