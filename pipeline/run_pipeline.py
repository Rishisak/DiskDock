
import json
from pathlib import Path
import re

from compound.structure import smiles_to_3d
from compound.preparation import prepare_ligand

from protein.structure import extract_chain
from protein.preparation import clean_protein, prepare_receptor

from docking.box import (
    calculate_box_from_reference_ligand,
    calculate_box_from_residues,
    calculate_box_from_coordinates,
    calculate_blind_docking_box
)

from docking.vina_runner import run_vina


def run_pipeline(
    protein_pdb,
    compound_smiles,
    protein_name="protein",
    compound_name="compound",
    chain_id="A",
    site_mode="blind",
    reference_ligand=None,
    residue_numbers=None,
    center=None,
    box_size=None,
    results_dir="results",
    vina_seed=42,
):
    """
    Run one complete protein-compound docking job.
    """

    def safe_name(value):
        cleaned = re.sub(
            r"[^A-Za-z0-9._-]+",
            "_",
            value.strip()
        )

        return cleaned.strip("._") or "unnamed"

    work_dir = (
        Path(results_dir)
        / safe_name(protein_name)
        / safe_name(compound_name)
    )

    work_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    print("\n========================================")
    print("       DistDockNet Docking Pipeline")
    print("========================================")

    # --------------------------------------------------
    # 1. Extract protein chain
    # --------------------------------------------------

    chain_pdb = work_dir / "protein_chain.pdb"

    chain_pdb, selected_chain_id = extract_chain(
        str(protein_pdb),
        str(chain_pdb),
        chain_id
    )

    # --------------------------------------------------
    # 2. Clean protein
    # --------------------------------------------------

    clean_pdb = work_dir / "protein_clean.pdb"

    clean_protein(
        str(chain_pdb),
        str(clean_pdb)
    )

    # --------------------------------------------------
    # 3. Prepare receptor
    # --------------------------------------------------

    receptor_prefix = work_dir / "receptor"

    prepare_receptor(
        str(clean_pdb),
        str(receptor_prefix)
    )

    receptor_pdbqt = str(
        receptor_prefix
    ) + ".pdbqt"

    # --------------------------------------------------
    # 4. Convert compound SMILES to 3D
    # --------------------------------------------------

    ligand_sdf = work_dir / "ligand.sdf"

    smiles_to_3d(
        compound_smiles,
        str(ligand_sdf)
    )

    # --------------------------------------------------
    # 5. Prepare ligand
    # --------------------------------------------------

    ligand_pdbqt = work_dir / "ligand.pdbqt"

    prepare_ligand(
        str(ligand_sdf),
        str(ligand_pdbqt)
    )

    # --------------------------------------------------
    # 6. Calculate docking box
    # --------------------------------------------------

    print("\n========================================")
    print("          DOCKING SITE")
    print("========================================")

    print(
        "Selected mode:",
        site_mode
    )

    if site_mode == "reference":

        if reference_ligand is None:
            raise ValueError(
                "Reference ligand was not provided."
            )

        box = calculate_box_from_reference_ligand(
            str(protein_pdb),
            reference_ligand,
            padding=10.0
        )

    elif site_mode == "residues":

        if not residue_numbers:
            raise ValueError(
                "Binding-site residues were not provided."
            )

        box = calculate_box_from_residues(
            str(protein_pdb),
            residue_numbers,
            chain_id=selected_chain_id,
            padding=8.0
        )

    elif site_mode == "coordinates":

        if center is None:
            raise ValueError(
                "Docking center was not provided."
            )

        if box_size is None:
            box_size = {
                "x": 20.0,
                "y": 20.0,
                "z": 20.0
            }

        box = calculate_box_from_coordinates(
            center_x=center["x"],
            center_y=center["y"],
            center_z=center["z"],
            size_x=box_size["x"],
            size_y=box_size["y"],
            size_z=box_size["z"]
        )

    elif site_mode == "blind":

        box = calculate_blind_docking_box(
            str(clean_pdb)
        )

    else:

        raise ValueError(
            f"Unknown docking site mode: {site_mode}"
        )

    # --------------------------------------------------
    # Force box values to normal Python floats
    # --------------------------------------------------

    box = {
        "center_x": float(box["center_x"]),
        "center_y": float(box["center_y"]),
        "center_z": float(box["center_z"]),
        "size_x": float(box["size_x"]),
        "size_y": float(box["size_y"]),
        "size_z": float(box["size_z"]),
    }

    # --------------------------------------------------
    # 7. Run AutoDock Vina
    # --------------------------------------------------

    docking_output = (
        work_dir / "docking_output.pdbqt"
    )

    vina_executable = (
        Path(__file__).resolve().parents[1]
        / "vina.exe"
    )

    vina_result = run_vina(
        str(vina_executable),
        str(receptor_pdbqt),
        str(ligand_pdbqt),
        str(docking_output),
        box,
        exhaustiveness=16,
        num_modes=10,
        seed=int(vina_seed),
    )

    # --------------------------------------------------
    # 8. Convert Vina result to normal Python types
    # --------------------------------------------------

    binding_affinity = float(
        vina_result["binding_affinity"]
    )

    # --------------------------------------------------
    # 9. Save metadata
    # --------------------------------------------------

    metadata = {
        "protein_name": str(protein_name),
        "compound_name": str(compound_name),
        "protein_pdb": str(
            Path(protein_pdb).resolve()
        ),
        "requested_chain_id": str(chain_id),
        "selected_chain_id": str(
            selected_chain_id
        ),
        "site_mode": str(site_mode),

        "box": {
            "center_x": float(
                box["center_x"]
            ),
            "center_y": float(
                box["center_y"]
            ),
            "center_z": float(
                box["center_z"]
            ),
            "size_x": float(
                box["size_x"]
            ),
            "size_y": float(
                box["size_y"]
            ),
            "size_z": float(
                box["size_z"]
            ),
        },

        "binding_affinity": float(
            binding_affinity
        ),

        "vina_seed": int(vina_seed),
    }

    metadata_file = (
        work_dir / "run_metadata.json"
    )

    metadata_file.write_text(
        json.dumps(
            metadata,
            indent=2
        ),
        encoding="utf-8"
    )

    # --------------------------------------------------
    # 10. Final success message
    # --------------------------------------------------

    print("\n========================================")
    print("       PIPELINE COMPLETED")
    print("========================================")

    print(
        "\nProtein :",
        protein_name
    )

    print(
        "Compound:",
        compound_name
    )

    print(
        "Affinity:",
        binding_affinity
    )

    print(
        "Output  :",
        docking_output
    )

    print(
        "Metadata:",
        metadata_file
    )

    return {
        "work_directory": str(
            work_dir
        ),

        "docking_output": str(
            docking_output
        ),

        "binding_affinity": float(
            binding_affinity
        ),

        "selected_chain_id": str(
            selected_chain_id
        ),
    }
