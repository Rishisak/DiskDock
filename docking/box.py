from Bio.PDB import PDBParser


def _make_box(
    center_x,
    center_y,
    center_z,
    size_x,
    size_y,
    size_z
):
    """
    Create a docking-box dictionary using normal Python floats.

    This is important because Biopython coordinates may be
    NumPy float32 values, which cannot always be written directly
    into JSON.
    """

    return {
        "center_x": float(center_x),
        "center_y": float(center_y),
        "center_z": float(center_z),
        "size_x": float(size_x),
        "size_y": float(size_y),
        "size_z": float(size_z),
    }


def get_ligand_coordinates(pdb_file, ligand_resname):
    """
    Find all atoms belonging to a reference ligand.
    """

    pdb_file = str(pdb_file)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(
        "structure",
        pdb_file
    )

    coordinates = []

    for model in structure:
        for chain in model:
            for residue in chain:

                if residue.resname.strip() == ligand_resname:

                    for atom in residue:
                        coordinates.append(atom.coord)

    return coordinates


def calculate_box_from_reference_ligand(
    pdb_file,
    ligand_resname,
    padding=10.0
):
    """
    Calculate docking box around a reference ligand.
    """

    print(
        "\n[11] Calculating docking box "
        "from reference ligand..."
    )

    coordinates = get_ligand_coordinates(
        pdb_file,
        ligand_resname
    )

    if not coordinates:
        raise ValueError(
            f"Reference ligand '{ligand_resname}' "
            "was not found in the protein structure."
        )

    xs = [float(coord[0]) for coord in coordinates]
    ys = [float(coord[1]) for coord in coordinates]
    zs = [float(coord[2]) for coord in coordinates]

    min_x = min(xs)
    max_x = max(xs)

    min_y = min(ys)
    max_y = max(ys)

    min_z = min(zs)
    max_z = max(zs)

    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    center_z = (min_z + max_z) / 2

    size_x = (max_x - min_x) + 2 * float(padding)
    size_y = (max_y - min_y) + 2 * float(padding)
    size_z = (max_z - min_z) + 2 * float(padding)

    box = _make_box(
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z
    )

    print("\nDocking box:")
    print(f"center_x = {box['center_x']:.3f}")
    print(f"center_y = {box['center_y']:.3f}")
    print(f"center_z = {box['center_z']:.3f}")

    print(f"size_x   = {box['size_x']:.3f}")
    print(f"size_y   = {box['size_y']:.3f}")
    print(f"size_z   = {box['size_z']:.3f}")

    return box


def calculate_box_from_residues(
    pdb_file,
    residue_numbers,
    chain_id="A",
    padding=8.0
):
    """
    Calculate docking box around specified protein residues.

    Example:
        residue_numbers = [198, 201, 202, 222]
    """

    print(
        "\n[11] Calculating docking box "
        "from binding-site residues..."
    )

    pdb_file = str(pdb_file)

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "protein",
        pdb_file
    )

    coordinates = []

    for model in structure:

        if chain_id not in model:
            continue

        chain = model[chain_id]

        for residue in chain:

            residue_number = residue.id[1]

            if residue_number in residue_numbers:

                for atom in residue:
                    coordinates.append(atom.coord)

    if not coordinates:
        raise ValueError(
            "None of the specified binding-site residues "
            "were found."
        )

    xs = [float(coord[0]) for coord in coordinates]
    ys = [float(coord[1]) for coord in coordinates]
    zs = [float(coord[2]) for coord in coordinates]

    min_x = min(xs)
    max_x = max(xs)

    min_y = min(ys)
    max_y = max(ys)

    min_z = min(zs)
    max_z = max(zs)

    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    center_z = (min_z + max_z) / 2

    size_x = (max_x - min_x) + 2 * float(padding)
    size_y = (max_y - min_y) + 2 * float(padding)
    size_z = (max_z - min_z) + 2 * float(padding)

    box = _make_box(
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z
    )

    print("\nDocking box:")
    print(f"center_x = {box['center_x']:.3f}")
    print(f"center_y = {box['center_y']:.3f}")
    print(f"center_z = {box['center_z']:.3f}")

    print(f"size_x   = {box['size_x']:.3f}")
    print(f"size_y   = {box['size_y']:.3f}")
    print(f"size_z   = {box['size_z']:.3f}")

    return box


def calculate_box_from_coordinates(
    center_x,
    center_y,
    center_z,
    size_x=20.0,
    size_y=20.0,
    size_z=20.0
):
    """
    Create docking box from user-provided coordinates.
    """

    print(
        "\n[11] Using user-provided "
        "docking coordinates..."
    )

    box = _make_box(
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z
    )

    print("\nDocking box:")
    print(f"center_x = {box['center_x']}")
    print(f"center_y = {box['center_y']}")
    print(f"center_z = {box['center_z']}")

    print(f"size_x   = {box['size_x']}")
    print(f"size_y   = {box['size_y']}")
    print(f"size_z   = {box['size_z']}")

    return box


def calculate_blind_docking_box(pdb_file):
    """
    Create a large docking box covering the protein.
    """

    print("\n[11] Calculating blind docking box...")

    pdb_file = str(pdb_file)

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure(
        "protein",
        pdb_file
    )

    coordinates = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.coord)

    if not coordinates:
        raise ValueError(
            "No protein coordinates found."
        )

    xs = [float(coord[0]) for coord in coordinates]
    ys = [float(coord[1]) for coord in coordinates]
    zs = [float(coord[2]) for coord in coordinates]

    min_x = min(xs)
    max_x = max(xs)

    min_y = min(ys)
    max_y = max(ys)

    min_z = min(zs)
    max_z = max(zs)

    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    center_z = (min_z + max_z) / 2

    padding = 5.0

    size_x = (max_x - min_x) + 2 * padding
    size_y = (max_y - min_y) + 2 * padding
    size_z = (max_z - min_z) + 2 * padding

    box = _make_box(
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z
    )

    print("\nBlind docking box:")
    print(f"center_x = {box['center_x']:.3f}")
    print(f"center_y = {box['center_y']:.3f}")
    print(f"center_z = {box['center_z']:.3f}")

    print(f"size_x   = {box['size_x']:.3f}")
    print(f"size_y   = {box['size_y']:.3f}")
    print(f"size_z   = {box['size_z']:.3f}")

    return box
