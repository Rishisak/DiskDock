from docking.box import calculate_box_from_reference_ligand


protein = "protein/6GU2.pdb"

box = calculate_box_from_reference_ligand(
    protein,
    "F9Z",
    padding=10.0
)

print("\nFinal box:")
print(box)