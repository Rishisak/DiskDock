import re
import subprocess
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
import streamlit as st
import py3Dmol
import streamlit.components.v1 as components

from rdkit import Chem
from rdkit.Chem import AllChem

from Bio.PDB import PDBParser, PDBIO, Select

from meeko import MoleculePreparation, PDBQTWriterLegacy

try:
    from analysis import (
        fetch_string_ppi,
        build_ppi_graph,
        calculate_network_centralities,
        normalize_target_to_gene,
        fetch_functional_enrichment,
        build_target_compound_network,
        rank_compounds_by_polypharmacology,
        generate_network_html,
        get_enrichment_altair_chart,
    )
    _ANALYSIS_OK = True
except ImportError:
    _ANALYSIS_OK = False


# ============================================================
# PAGE CONFIG
# ============================================================

st.set_page_config(
    page_title="DistDockNet",
    page_icon="🧬",
    layout="wide",
)


# ============================================================
# PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parent

VINA_EXE = BASE_DIR / "vina.exe"

RESULTS_DIR = (
    BASE_DIR
    / "results"
    / "ui_runs"
)

PROTEIN_CACHE = (
    BASE_DIR
    / "protein"
    / "cache"
)

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

PROTEIN_CACHE.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# STANDARD AMINO ACIDS
# ============================================================

STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


# ============================================================
# GENERAL HELPERS
# ============================================================

def safe_name(value):
    value = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(value).strip(),
    )

    value = value[:80]

    return value or "item"


# ============================================================
# RCSB PDB DOWNLOAD
# ============================================================

def fetch_pdb(pdb_id):
    import urllib.request

    pdb_id = str(pdb_id).strip().upper()

    if not re.fullmatch(
        r"[0-9A-Za-z]{4}",
        pdb_id,
    ):
        raise ValueError(
            "PDB ID must contain exactly 4 letters/numbers."
        )

    output = (
        PROTEIN_CACHE
        / f"{pdb_id}.pdb"
    )

    # Use cached copy if available.
    if (
        output.exists()
        and output.stat().st_size > 0
    ):
        return output

    url = (
        f"https://files.rcsb.org/download/"
        f"{pdb_id}.pdb"
    )

    try:
        urllib.request.urlretrieve(
            url,
            output,
        )

    except Exception as exc:

        if output.exists():
            output.unlink()

        raise RuntimeError(
            f"Could not download PDB {pdb_id}: {exc}"
        ) from exc

    if (
        not output.exists()
        or output.stat().st_size == 0
    ):
        raise RuntimeError(
            f"RCSB returned an empty PDB file for {pdb_id}."
        )

    return output


# ============================================================
# SMILES -> 3D SDF
# ============================================================

def smiles_to_sdf(
    smiles,
    output_sdf,
):
    """
    Convert SMILES into an optimized 3D SDF file.
    """

    smiles = str(smiles).strip()

    if not smiles:
        raise ValueError(
            "SMILES cannot be empty."
        )

    mol = Chem.MolFromSmiles(
        smiles
    )

    if mol is None:
        raise ValueError(
            "Invalid SMILES string."
        )

    mol = Chem.AddHs(
        mol
    )

    params = AllChem.ETKDGv3()

    params.randomSeed = 42

    status = AllChem.EmbedMolecule(
        mol,
        params,
    )

    if status != 0:

        # Fallback.
        status = AllChem.EmbedMolecule(
            mol,
            randomSeed=42,
        )

    if status != 0:
        raise RuntimeError(
            "RDKit could not generate "
            "a 3D conformer."
        )

    try:

        AllChem.UFFOptimizeMolecule(
            mol
        )

    except Exception:

        # The molecule can still be used
        # if UFF optimization fails.
        pass

    output_sdf = Path(
        output_sdf
    )

    output_sdf.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    writer = Chem.SDWriter(
        str(output_sdf)
    )

    writer.write(
        mol
    )

    writer.close()

    if (
        not output_sdf.exists()
        or output_sdf.stat().st_size == 0
    ):
        raise RuntimeError(
            "RDKit did not create the SDF file."
        )

    return output_sdf


# ============================================================
# PROTEIN SELECTOR
# ============================================================

class ProteinSelect(Select):

    def __init__(
        self,
        chain_id,
    ):
        self.chain_id = chain_id

    def accept_chain(
        self,
        chain,
    ):
        return (
            chain.id
            == self.chain_id
        )

    def accept_residue(
        self,
        residue,
    ):

        if (
            residue.resname.strip()
            not in STANDARD_AMINO_ACIDS
        ):
            return 0

        if residue.id[1] < 0:
            return 0

        return 1


# ============================================================
# EXTRACT + CLEAN PROTEIN
# ============================================================

def extract_and_clean_protein(
    input_pdb,
    output_pdb,
    chain_id="A",
):

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "protein",
        str(input_pdb),
    )

    model = next(
        structure.get_models(),
        None,
    )

    if model is None:
        raise ValueError(
            "Protein has no model."
        )

    available = {}

    for chain in model:

        count = sum(
            1
            for residue in chain
            if residue.resname.strip()
            in STANDARD_AMINO_ACIDS
        )

        if count:
            available[
                chain.id
            ] = count

    if not available:
        raise ValueError(
            "No protein chains found."
        )

    if chain_id in available:

        selected = chain_id

    else:

        selected = max(
            available,
            key=available.get,
        )

    output_pdb = Path(
        output_pdb
    )

    output_pdb.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    io_obj = PDBIO()

    io_obj.set_structure(
        structure
    )

    io_obj.save(
        str(output_pdb),
        ProteinSelect(
            selected
        ),
    )

    if (
        not output_pdb.exists()
        or output_pdb.stat().st_size == 0
    ):
        raise RuntimeError(
            "Clean protein PDB was not created."
        )

    return (
        selected,
        available,
    )


# ============================================================
# EXTRACT EXPERIMENTAL LIGAND
# ============================================================

def extract_reference_ligand(
    input_pdb,
    output_pdb,
    ligand_name,
    chain_id=None,
):

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "complex",
        str(input_pdb),
    )

    matches = []

    for model in structure:

        for chain in model:

            if (
                chain_id
                and chain.id != chain_id
            ):
                continue

            for residue in chain:

                if (
                    residue.resname.strip()
                    == ligand_name.strip()
                ):

                    matches.append(
                        (
                            model.id,
                            chain.id,
                            residue.id,
                        )
                    )

    if not matches:

        raise ValueError(
            f"Reference ligand "
            f"'{ligand_name}' "
            f"was not found in the PDB."
        )

    model_id, selected_chain, selected_residue_id = (
        matches[0]
    )

    class SingleLigandSelect(
        Select
    ):

        def accept_model(
            self,
            model,
        ):
            return (
                model.id
                == model_id
            )

        def accept_chain(
            self,
            chain,
        ):
            return (
                chain.id
                == selected_chain
            )

        def accept_residue(
            self,
            residue,
        ):
            return (
                residue.resname.strip()
                == ligand_name.strip()
                and residue.id
                == selected_residue_id
            )

    output_pdb = Path(
        output_pdb
    )

    output_pdb.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    io_obj = PDBIO()

    io_obj.set_structure(
        structure
    )

    io_obj.save(
        str(output_pdb),
        SingleLigandSelect(),
    )

    return selected_chain


# ============================================================
# LOAD REFERENCE LIGAND
# ============================================================

def load_reference_ligand_as_rdkit(
    ligand_pdb,
):

    mol = Chem.MolFromPDBFile(
        str(ligand_pdb),
        removeHs=False,
        sanitize=False,
    )

    if mol is None:
        raise ValueError(
            "RDKit could not read "
            "the extracted reference ligand."
        )

    try:

        Chem.SanitizeMol(
            mol
        )

    except Exception:

        pass

    if (
        mol.GetNumConformers()
        == 0
    ):
        raise ValueError(
            "Reference ligand has "
            "no 3D coordinates."
        )

    return mol


# ============================================================
# MEeko LIGAND PREPARATION
# ============================================================

def prepare_rdkit_ligand(
    mol,
    output_pdbqt,
):

    if mol is None:
        raise ValueError(
            "Cannot prepare a None ligand."
        )

    mol = Chem.Mol(
        mol
    )

    try:

        mol = Chem.AddHs(
            mol,
            addCoords=True,
        )

    except Exception as exc:

        raise RuntimeError(
            "Could not add explicit "
            f"hydrogens: {exc}"
        ) from exc

    if (
        mol.GetNumConformers()
        == 0
    ):

        params = AllChem.ETKDGv3()

        params.randomSeed = 42

        status = AllChem.EmbedMolecule(
            mol,
            params,
        )

        if status != 0:

            raise RuntimeError(
                "Could not generate "
                "ligand 3D coordinates."
            )

    preparator = MoleculePreparation()

    try:

        setups = preparator.prepare(
            mol
        )

    except Exception as exc:

        raise RuntimeError(
            "Meeko ligand preparation failed: "
            f"{exc}"
        ) from exc

    if not setups:
        raise RuntimeError(
            "Meeko did not generate "
            "a ligand setup."
        )

    pdbqt_string, is_ok, error_msg = (
        PDBQTWriterLegacy.write_string(
            setups[0]
        )
    )

    if not is_ok:

        raise RuntimeError(
            "Meeko PDBQT writing failed: "
            f"{error_msg}"
        )

    output_pdbqt = Path(
        output_pdbqt
    )

    output_pdbqt.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_pdbqt.write_text(
        pdbqt_string,
        encoding="utf-8",
    )

    if (
        not output_pdbqt.exists()
        or output_pdbqt.stat().st_size == 0
    ):

        raise RuntimeError(
            "Ligand PDBQT was not created."
        )

    return output_pdbqt


def prepare_ligand(
    sdf_file,
    output_pdbqt,
):

    supplier = Chem.SDMolSupplier(
        str(sdf_file),
        removeHs=False,
    )

    if (
        not supplier
        or supplier[0] is None
    ):

        raise ValueError(
            "Could not read "
            "generated ligand structure."
        )

    return prepare_rdkit_ligand(
        supplier[0],
        output_pdbqt,
    )


# ============================================================
# MEeko RECEPTOR PREPARATION
# ============================================================

def prepare_receptor(
    pdb_file,
    output_prefix,
):

    command = [
        "mk_prepare_receptor.exe",
        "--read_pdb",
        str(pdb_file),
        "-o",
        str(output_prefix),
        "-p",
        "--default_altloc",
        "A",
        "--allow_bad_res",
    ]

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

    except OSError as exc:

        raise RuntimeError(
            "Could not start "
            "mk_prepare_receptor.exe. "
            "Make sure Meeko is installed "
            "and the executable is available."
        ) from exc

    if result.returncode != 0:

        details = (
            result.stderr
            or result.stdout
            or "Unknown Meeko error."
        )

        raise RuntimeError(
            "Meeko receptor preparation failed.\n\n"
            + details
        )

    receptor = Path(
        str(output_prefix)
        + ".pdbqt"
    )

    if not receptor.exists():

        raise RuntimeError(
            "Meeko completed but receptor "
            "PDBQT was not created."
        )

    return receptor


# ============================================================
# PROTEIN COORDINATES
# ============================================================

def protein_coordinates(
    pdb_file,
):

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "protein",
        str(pdb_file),
    )

    model = next(
        structure.get_models(),
        None,
    )

    if model is None:
        raise ValueError(
            "Protein has no model."
        )

    coordinates = []

    for chain in model:

        for residue in chain:

            if (
                residue.resname.strip()
                not in STANDARD_AMINO_ACIDS
            ):
                continue

            for atom in residue:

                coordinates.append(
                    atom.coord
                )

    if not coordinates:
        raise ValueError(
            "No protein coordinates found."
        )

    return coordinates


# ============================================================
# BLIND DOCKING BOX
# ============================================================

def calculate_blind_box(
    pdb_file,
):

    coordinates = protein_coordinates(
        pdb_file
    )

    arr = np.asarray(
        coordinates,
        dtype=float,
    )

    minimum = arr.min(
        axis=0
    )

    maximum = arr.max(
        axis=0
    )

    padding = 5.0

    size = (
        maximum
        - minimum
        + 2 * padding
    )

    center = (
        minimum
        + maximum
    ) / 2

    return {
        "center_x": float(
            center[0]
        ),
        "center_y": float(
            center[1]
        ),
        "center_z": float(
            center[2]
        ),
        "size_x": float(
            size[0]
        ),
        "size_y": float(
            size[1]
        ),
        "size_z": float(
            size[2]
        ),
    }


# ============================================================
# REFERENCE LIGAND DOCKING BOX
# ============================================================

def calculate_reference_ligand_box(
    pdb_file,
    ligand_resname,
    padding=8.0,
    chain_id=None,
):

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "protein",
        str(pdb_file),
    )

    model = next(
        structure.get_models(),
        None,
    )

    if model is None:
        raise ValueError(
            "Protein has no model."
        )

    coordinates = []

    for chain in model:

        if (
            chain_id
            and chain.id != chain_id
        ):
            continue

        for residue in chain:

            if (
                residue.resname.strip()
                == ligand_resname.strip()
            ):

                for atom in residue:

                    coordinates.append(
                        atom.coord
                    )

    if not coordinates:

        raise ValueError(
            f"Reference ligand "
            f"'{ligand_resname}' "
            "was not found."
        )

    arr = np.asarray(
        coordinates,
        dtype=float,
    )

    minimum = arr.min(
        axis=0
    )

    maximum = arr.max(
        axis=0
    )

    size = (
        maximum
        - minimum
        + 2 * padding
    )

    center = (
        minimum
        + maximum
    ) / 2

    return {
        "center_x": float(
            center[0]
        ),
        "center_y": float(
            center[1]
        ),
        "center_z": float(
            center[2]
        ),
        "size_x": float(
            size[0]
        ),
        "size_y": float(
            size[1]
        ),
        "size_z": float(
            size[2]
        ),
    }


# ============================================================
# RUN AUTODOCK VINA
# ============================================================

def run_vina(
    receptor,
    ligand,
    box,
    output_file,
    exhaustiveness=16,
    num_modes=10,
):

    output_file = Path(
        output_file
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    config_file = (
        output_file.parent
        / "vina_config.txt"
    )

    config_text = f"""
receptor = {Path(receptor).resolve()}
ligand = {Path(ligand).resolve()}

center_x = {box["center_x"]}
center_y = {box["center_y"]}
center_z = {box["center_z"]}

size_x = {box["size_x"]}
size_y = {box["size_y"]}
size_z = {box["size_z"]}

exhaustiveness = {int(exhaustiveness)}
num_modes = {int(num_modes)}

seed = 42
"""

    config_file.write_text(
        config_text.strip()
        + "\n",
        encoding="utf-8",
    )

    command = [
        str(VINA_EXE),
        "--config",
        str(config_file),
        "--out",
        str(output_file),
    ]

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

    except OSError as exc:

        raise RuntimeError(
            f"Could not start Vina: {exc}"
        ) from exc

    if result.returncode != 0:

        details = (
            result.stderr
            or result.stdout
            or "Unknown Vina error."
        )

        raise RuntimeError(
            "Vina docking failed.\n\n"
            + details
        )

    if (
        not output_file.exists()
        or output_file.stat().st_size == 0
    ):

        raise RuntimeError(
            "Vina finished but did not "
            "create a docking output file."
        )

    return result.stdout


# ============================================================
# PARSE VINA SCORE
# ============================================================

def parse_vina_score(
    pdbqt_file,
):

    pattern = re.compile(
        r"REMARK VINA RESULT:\s+"
        r"(-?\d+(?:\.\d+)?)"
    )

    text = Path(
        pdbqt_file
    ).read_text(
        encoding="utf-8",
        errors="ignore",
    )

    for line in text.splitlines():

        match = pattern.search(
            line
        )

        if match:

            return float(
                match.group(1)
            )

    raise ValueError(
        "Could not find VINA RESULT "
        f"in docking output: {pdbqt_file}"
    )


# ============================================================
# PDBQT ELEMENT DETECTION
# ============================================================

def _pdbqt_element(
    line,
):

    atom_type = ""

    if len(line) >= 79:

        atom_type = (
            line[77:79]
            .strip()
            .upper()
        )

    mapping = {
        "C": "C",
        "A": "C",
        "N": "N",
        "NA": "N",
        "NS": "N",
        "OA": "O",
        "O": "O",
        "OS": "O",
        "S": "S",
        "SA": "S",
        "P": "P",
        "F": "F",
        "CL": "CL",
        "BR": "BR",
        "I": "I",
        "H": "H",
        "HD": "H",
    }

    if atom_type in mapping:

        return mapping[
            atom_type
        ]

    atom_name = (
        line[12:16]
        .strip()
        .upper()
    )

    if atom_name.startswith(
        ("CL", "BR")
    ):

        return atom_name[:2]

    cleaned = re.sub(
        r"[^A-Z]",
        "",
        atom_name,
    )

    return (
        cleaned[:1]
        if cleaned
        else ""
    )


# ============================================================
# PARSE PDBQT ATOMS
# ============================================================

def parse_pdbqt_atoms(
    pdbqt_file,
    first_model_only=True,
):

    text = Path(
        pdbqt_file
    ).read_text(
        encoding="utf-8",
        errors="ignore",
    )

    lines = text.splitlines()

    atoms = []

    has_models = any(
        line.startswith("MODEL")
        for line in lines
    )

    in_model = not has_models

    model_seen = False

    for line in lines:

        if line.startswith(
            "MODEL"
        ):

            model_number = (
                line[5:]
                .strip()
            )

            if first_model_only:

                if model_seen:
                    break

                if model_number not in (
                    "1",
                    "01",
                ):
                    continue

                in_model = True
                model_seen = True

            else:

                in_model = True

        elif line.startswith(
            "ENDMDL"
        ):

            if (
                first_model_only
                and model_seen
            ):

                break

            in_model = False

        elif (
            in_model
            and line.startswith(
                (
                    "ATOM",
                    "HETATM",
                )
            )
        ):

            try:

                atom_name = (
                    line[12:16]
                    .strip()
                )

                x = float(
                    line[30:38]
                )

                y = float(
                    line[38:46]
                )

                z = float(
                    line[46:54]
                )

                atoms.append(
                    {
                        "name": atom_name,
                        "element": _pdbqt_element(
                            line
                        ),
                        "coord": np.array(
                            [
                                x,
                                y,
                                z,
                            ],
                            dtype=float,
                        ),
                    }
                )

            except (
                ValueError,
                IndexError,
            ):

                continue

    return atoms


# ============================================================
# PDB LIGAND ATOMS
# ============================================================

def parse_pdb_ligand_atoms_from_file(
    ligand_pdb,
):

    parser = PDBParser(
        QUIET=True
    )

    structure = parser.get_structure(
        "ligand",
        str(ligand_pdb),
    )

    model = next(
        structure.get_models(),
        None,
    )

    if model is None:
        return []

    atoms = []

    for chain in model:

        for residue in chain:

            if (
                residue.resname.strip()
                == "HOH"
            ):
                continue

            for atom in residue:

                element = (
                    atom.element
                    or ""
                ).strip().upper()

                if not element:

                    element = re.sub(
                        r"[^A-Za-z]",
                        "",
                        atom.name,
                    ).upper()[:1]

                atoms.append(
                    {
                        "name": atom.name.strip(),
                        "element": element,
                        "coord": np.array(
                            atom.coord,
                            dtype=float,
                        ),
                    }
                )

    return atoms


# ============================================================
# KABSCH RMSD
# ============================================================

def kabsch_rmsd(
    reference_coords,
    predicted_coords,
):

    reference = np.asarray(
        reference_coords,
        dtype=float,
    )

    predicted = np.asarray(
        predicted_coords,
        dtype=float,
    )

    if (
        reference.shape
        != predicted.shape
        or len(reference) < 3
    ):

        raise ValueError(
            "At least 3 matching "
            "atoms are required "
            "for RMSD."
        )

    ref_center = (
        reference.mean(
            axis=0
        )
    )

    pred_center = (
        predicted.mean(
            axis=0
        )
    )

    ref = (
        reference
        - ref_center
    )

    pred = (
        predicted
        - pred_center
    )

    covariance = (
        pred.T @ ref
    )

    v, _, wt = np.linalg.svd(
        covariance
    )

    determinant = np.linalg.det(
        v @ wt
    )

    d = np.sign(
        determinant
    )

    rotation = (
        v
        @ np.diag(
            [
                1.0,
                1.0,
                d,
            ]
        )
        @ wt
    )

    aligned = (
        pred
        @ rotation
    )

    difference = (
        ref
        - aligned
    )

    return float(
        np.sqrt(
            np.mean(
                np.sum(
                    difference ** 2,
                    axis=1,
                )
            )
        )
    )


# ============================================================
# REDOCKING RMSD
# ============================================================

def calculate_redocking_rmsd(
    reference_ligand_pdb,
    docked_pdbqt,
):

    """
    Compare experimental ligand coordinates
    with Vina pose 1.

    Only heavy atoms are compared.
    """

    reference_atoms = (
        parse_pdb_ligand_atoms_from_file(
            reference_ligand_pdb
        )
    )

    docked_atoms = (
        parse_pdbqt_atoms(
            docked_pdbqt,
            first_model_only=True,
        )
    )

    # Remove hydrogens.
    reference_atoms = [
        atom
        for atom in reference_atoms
        if atom["element"] != "H"
    ]

    docked_atoms = [
        atom
        for atom in docked_atoms
        if atom["element"] != "H"
    ]

    if (
        len(reference_atoms) < 3
        or len(docked_atoms) < 3
    ):

        raise ValueError(
            "Not enough heavy atoms "
            "available for RMSD."
        )

    reference_elements = [
        atom["element"]
        for atom in reference_atoms
    ]

    docked_elements = [
        atom["element"]
        for atom in docked_atoms
    ]

    pairs = []

    # --------------------------------------------------------
    # Best case:
    # same heavy atom count and same element order.
    # --------------------------------------------------------

    if (
        len(reference_atoms)
        == len(docked_atoms)
        and reference_elements
        == docked_elements
    ):

        pairs = list(
            zip(
                range(
                    len(reference_atoms)
                ),
                range(
                    len(docked_atoms)
                ),
            )
        )

    # --------------------------------------------------------
    # Second case:
    # match using atom names + element.
    # --------------------------------------------------------

    else:

        dock_by_name = {}

        for index, atom in enumerate(
            docked_atoms
        ):

            dock_by_name.setdefault(
                atom["name"].upper(),
                [],
            ).append(
                index
            )

        used = set()

        for ref_index, atom in enumerate(
            reference_atoms
        ):

            candidates = (
                dock_by_name.get(
                    atom["name"].upper(),
                    [],
                )
            )

            candidate = None

            for dock_index in candidates:

                if dock_index in used:
                    continue

                if (
                    docked_atoms[
                        dock_index
                    ]["element"]
                    == atom["element"]
                ):

                    candidate = (
                        dock_index
                    )

                    break

            if candidate is not None:

                pairs.append(
                    (
                        ref_index,
                        candidate,
                    )
                )

                used.add(
                    candidate
                )

    # --------------------------------------------------------
    # Third fallback:
    # match by element order.
    # --------------------------------------------------------

    if len(pairs) < 3:

        pairs = []

        element_positions = {}

        for index, element in enumerate(
            docked_elements
        ):

            element_positions.setdefault(
                element,
                [],
            ).append(
                index
            )

        counters = {}

        for ref_index, element in enumerate(
            reference_elements
        ):

            position = counters.get(
                element,
                0,
            )

            choices = (
                element_positions.get(
                    element,
                    [],
                )
            )

            if position < len(
                choices
            ):

                pairs.append(
                    (
                        ref_index,
                        choices[position],
                    )
                )

                counters[
                    element
                ] = (
                    position + 1
                )

    if len(pairs) < 3:

        raise ValueError(
            "Could not match enough "
            "heavy atoms for RMSD. "
            f"Reference heavy atoms: "
            f"{len(reference_atoms)}, "
            f"Docked heavy atoms: "
            f"{len(docked_atoms)}."
        )

    reference_coords = np.array(
        [
            reference_atoms[
                i
            ]["coord"]
            for i, _ in pairs
        ],
        dtype=float,
    )

    docked_coords = np.array(
        [
            docked_atoms[
                j
            ]["coord"]
            for _, j in pairs
        ],
        dtype=float,
    )

    rmsd = kabsch_rmsd(
        reference_coords,
        docked_coords,
    )

    return (
        float(rmsd),
        len(pairs),
        "heavy-atom correspondence + Kabsch",
    )


# ============================================================
# STANDARD DOCKING
# ============================================================

def run_standard_docking(
    protein_name,
    pdb_file,
    compound_name,
    smiles,
    chain_id,
    site_mode,
    reference_ligand,
    exhaustiveness,
    num_modes,
):

    run_dir = (
        RESULTS_DIR
        / safe_name(
            protein_name
        )
        / safe_name(
            compound_name
        )
    )

    run_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    clean_pdb = (
        run_dir
        / "protein_clean.pdb"
    )

    sdf_file = (
        run_dir
        / "ligand_3d.sdf"
    )

    ligand_pdbqt = (
        run_dir
        / "ligand.pdbqt"
    )

    receptor_prefix = (
        run_dir
        / "receptor"
    )

    output_pdbqt = (
        run_dir
        / "docking_output.pdbqt"
    )

    selected_chain, _ = (
        extract_and_clean_protein(
            pdb_file,
            clean_pdb,
            chain_id,
        )
    )

    smiles_to_sdf(
        smiles,
        sdf_file,
    )

    prepare_ligand(
        sdf_file,
        ligand_pdbqt,
    )

    receptor = (
        prepare_receptor(
            clean_pdb,
            receptor_prefix,
        )
    )

    if (
        site_mode
        == "reference"
    ):

        box = (
            calculate_reference_ligand_box(
                pdb_file,
                reference_ligand,
                padding=8.0,
                chain_id=selected_chain,
            )
        )

    else:

        box = calculate_blind_box(
            clean_pdb
        )

    vina_output = run_vina(
        receptor,
        ligand_pdbqt,
        box,
        output_pdbqt,
        exhaustiveness,
        num_modes,
    )

    score = parse_vina_score(
        output_pdbqt
    )

    return {
        "protein": protein_name,
        "compound": compound_name,
        "binding_affinity": score,
        "selected_chain": selected_chain,
        "output": str(
            output_pdbqt
        ),
        "protein_pdb": str(
            clean_pdb
        ),
        "ligand_pdbqt": str(
            ligand_pdbqt
        ),
        "box": box,
        "vina_output": vina_output,
    }


# ============================================================
# REDOCKING VALIDATION
# ============================================================

def run_redocking_validation(
    protein_name,
    pdb_file,
    reference_ligand,
    chain_id,
    exhaustiveness,
    num_modes,
):

    run_dir = (
        RESULTS_DIR
        / safe_name(
            protein_name
        )
        / (
            "redocking_"
            + safe_name(
                reference_ligand
            )
        )
    )

    run_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    clean_pdb = (
        run_dir
        / "protein_clean.pdb"
    )

    reference_ligand_pdb = (
        run_dir
        / "reference_ligand_experimental.pdb"
    )

    reference_ligand_pdbqt = (
        run_dir
        / "reference_ligand.pdbqt"
    )

    receptor_prefix = (
        run_dir
        / "receptor"
    )

    output_pdbqt = (
        run_dir
        / "redocking_output.pdbqt"
    )

    selected_chain, _ = (
        extract_and_clean_protein(
            pdb_file,
            clean_pdb,
            chain_id,
        )
    )

    extract_reference_ligand(
        pdb_file,
        reference_ligand_pdb,
        reference_ligand,
        selected_chain,
    )

    reference_mol = (
        load_reference_ligand_as_rdkit(
            reference_ligand_pdb
        )
    )

    prepare_rdkit_ligand(
        reference_mol,
        reference_ligand_pdbqt,
    )

    receptor = (
        prepare_receptor(
            clean_pdb,
            receptor_prefix,
        )
    )

    box = (
        calculate_reference_ligand_box(
            pdb_file,
            reference_ligand,
            padding=8.0,
            chain_id=selected_chain,
        )
    )

    vina_output = run_vina(
        receptor,
        reference_ligand_pdbqt,
        box,
        output_pdbqt,
        exhaustiveness,
        num_modes,
    )

    score = parse_vina_score(
        output_pdbqt
    )

    # IMPORTANT:
    # calculate_redocking_rmsd accepts exactly
    # two arguments.
    rmsd, matched_atoms, rmsd_method = (
        calculate_redocking_rmsd(
            reference_ligand_pdb,
            output_pdbqt,
        )
    )

    return {
        "protein": protein_name,
        "reference_ligand": reference_ligand,
        "binding_affinity": score,
        "rmsd": rmsd,
        "matched_atoms": matched_atoms,
        "rmsd_method": rmsd_method,
        "selected_chain": selected_chain,
        "experimental_ligand": str(
            reference_ligand_pdb
        ),
        "predicted_output": str(
            output_pdbqt
        ),
        "protein_pdb": str(
            clean_pdb
        ),
        "receptor_pdbqt": str(
            receptor
        ),
        "reference_ligand_pdbqt": str(
            reference_ligand_pdbqt
        ),
        "box": box,
        "vina_output": vina_output,
    }


# ============================================================
# 3D VISUALIZATION
# ============================================================

def _read_text_file(
    path,
):

    return Path(
        path
    ).read_text(
        encoding="utf-8",
        errors="ignore",
    )


# ============================================================
# EXTRACT VINA POSES
# ============================================================

def _extract_vina_poses(
    pdbqt_text,
):

    poses = []

    model_matches = re.findall(
        r"MODEL\s+(\d+)"
        r"(.*?)"
        r"(?:ENDMDL|$)",
        pdbqt_text,
        re.DOTALL,
    )

    for (
        model_number,
        model_text,
    ) in model_matches:

        score_match = re.search(
            r"REMARK VINA RESULT:\s+"
            r"(-?\d+(?:\.\d+)?)\s+"
            r"(-?\d+(?:\.\d+)?)\s+"
            r"(-?\d+(?:\.\d+)?)",
            model_text,
        )

        if score_match:

            score = float(
                score_match.group(1)
            )

            rmsd_lb = float(
                score_match.group(2)
            )

            rmsd_ub = float(
                score_match.group(3)
            )

        else:

            score = None
            rmsd_lb = None
            rmsd_ub = None

        poses.append(
            {
                "pose": int(
                    model_number
                ),
                "score": score,
                "rmsd_lb": rmsd_lb,
                "rmsd_ub": rmsd_ub,
                "pdbqt": (
                    "MODEL\n"
                    + model_text
                    + "\nENDMDL\n"
                ),
            }
        )

    # Fallback for single-pose output.
    if (
        not poses
        and pdbqt_text.strip()
    ):

        score_match = re.search(
            r"REMARK VINA RESULT:\s+"
            r"(-?\d+(?:\.\d+)?)",
            pdbqt_text,
        )

        poses.append(
            {
                "pose": 1,
                "score": (
                    float(
                        score_match.group(1)
                    )
                    if score_match
                    else None
                ),
                "rmsd_lb": 0.0,
                "rmsd_ub": 0.0,
                "pdbqt": pdbqt_text,
            }
        )

    return poses


# ============================================================
# CONVERT SELECTED PDBQT POSE -> PDB
# ============================================================

def pdbqt_pose_to_pdb(
    pdbqt_text,
):

    lines = []

    serial = 1

    for line in pdbqt_text.splitlines():

        if not line.startswith(
            (
                "ATOM",
                "HETATM",
            )
        ):
            continue

        try:

            atom_name = (
                line[12:16]
                .strip()
                or "X"
            )

            x = float(
                line[30:38]
            )

            y = float(
                line[38:46]
            )

            z = float(
                line[46:54]
            )

            element = _pdbqt_element(
                line
            )

            if element == "CL":
                display_element = "Cl"

            elif element == "BR":
                display_element = "Br"

            else:
                display_element = (
                    element[:1]
                    .upper()
                )

            pdb_line = (
                f"HETATM"
                f"{serial:5d} "
                f"{atom_name:<4s}"
                f" LIG A   1"
                f"{x:12.3f}"
                f"{y:8.3f}"
                f"{z:8.3f}"
                f"{1.00:6.2f}"
                f"{0.00:6.2f}"
                f"          "
                f"{display_element:>2s}"
            )

            lines.append(
                pdb_line
            )

            serial += 1

        except (
            ValueError,
            IndexError,
        ):

            continue

    if not lines:
        return ""

    return (
        "\n".join(
            lines
        )
        + "\nEND\n"
    )


# ============================================================
# CREATE 3D VIEWER
# ============================================================

def _create_docking_viewer(
    protein_data,
    ligand_data,
    protein_style="Cartoon",
    ligand_style="Sticks",
    dark_background=True,
):

    viewer = py3Dmol.view(
        width="100%",
        height=680,
    )

    # --------------------------------------------------------
    # PROTEIN
    # --------------------------------------------------------

    viewer.addModel(
        protein_data,
        "pdb",
    )

    if protein_style == "Cartoon":

        viewer.setStyle(
            {"model": 0},
            {
                "cartoon": {
                    "color": "spectrum",
                    "opacity": 0.88,
                }
            },
        )

    elif protein_style == "Surface":

        viewer.setStyle(
            {"model": 0},
            {
                "cartoon": {
                    "color": "spectrum",
                    "opacity": 0.55,
                }
            },
        )

        viewer.addSurface(
            py3Dmol.VDW,
            {
                "opacity": 0.24,
                "color": "white",
            },
            {"model": 0},
        )

    else:

        viewer.setStyle(
            {"model": 0},
            {
                "stick": {
                    "radius": 0.12,
                }
            },
        )

    # --------------------------------------------------------
    # LIGAND
    # --------------------------------------------------------

    viewer.addModel(
        ligand_data,
        "pdb",
    )

    if ligand_style == "Sticks":

        viewer.setStyle(
            {"model": 1},
            {
                "stick": {
                    "radius": 0.25,
                }
            },
        )

    elif ligand_style == "Ball & Stick":

        viewer.setStyle(
            {"model": 1},
            {
                "stick": {
                    "radius": 0.20,
                },
                "sphere": {
                    "scale": 0.30,
                },
            },
        )

    else:

        viewer.setStyle(
            {"model": 1},
            {
                "sphere": {
                    "scale": 0.45,
                }
            },
        )

    # --------------------------------------------------------
    # VIEW SETTINGS
    # --------------------------------------------------------

    if dark_background:

        viewer.setBackgroundColor(
            "#0b1020"
        )

    else:

        viewer.setBackgroundColor(
            "white"
        )

    viewer.setProjection(
        "perspective"
    )

    viewer.zoomTo()

    viewer.zoom(
        1.15
    )

    viewer.spin(
        False
    )

    return viewer


# ============================================================
# RENDER 3D VIEWER
# ============================================================

def _render_docking_3d_viewer(
    results_df,
):

    if (
        results_df is None
        or results_df.empty
    ):
        return

    affinity_column = (
        "Binding affinity (kcal/mol)"
    )

    if (
        affinity_column
        not in results_df.columns
    ):
        return

    numeric_affinity = pd.to_numeric(
        results_df[
            affinity_column
        ],
        errors="coerce",
    )

    successful = results_df.loc[
        numeric_affinity.notna()
    ].copy()

    if successful.empty:
        return

    st.divider()

    st.subheader(
        "5. Interactive 3D docking visualization"
    )

    st.caption(
        "Select a completed docking job "
        "and inspect the protein together "
        "with one of the poses produced "
        "by AutoDock Vina."
    )

    # --------------------------------------------------------
    # SELECT DOCKING JOB
    # --------------------------------------------------------

    options = list(
        successful.index
    )

    labels = {}

    for index in options:

        row = successful.loc[
            index
        ]

        labels[index] = (
            f"{row['Protein']} + "
            f"{row['Compound']} "
            f"("
            f"{float(row[affinity_column]):.3f}"
            f" kcal/mol)"
        )

    selected_index = st.selectbox(
        "Docking complex",
        options,
        format_func=lambda index:
            labels[index],
        key="three_d_result_selection",
    )

    selected = successful.loc[
        selected_index
    ]

    # --------------------------------------------------------
    # DOCKING OUTPUT
    # --------------------------------------------------------

    output_value = str(
        selected.get(
            "Docking output",
            "",
        )
        or ""
    ).strip()

    if not output_value:

        st.error(
            "This docking result does not "
            "contain a docking output path."
        )

        return

    output_path = Path(
        output_value
    )

    if not output_path.is_file():

        st.error(
            "Docking output file was not found."
        )

        st.code(
            str(output_path)
        )

        st.info(
            "Run the docking again so "
            "that docking_output.pdbqt "
            "is generated."
        )

        return

    # --------------------------------------------------------
    # PROTEIN
    # --------------------------------------------------------

    protein_value = str(
        selected.get(
            "Protein PDB",
            "",
        )
        or ""
    ).strip()

    if protein_value:

        protein_path = Path(
            protein_value
        )

    else:

        protein_path = (
            output_path.parent
            / "protein_clean.pdb"
        )

    if not protein_path.is_file():

        fallback = (
            output_path.parent
            / "protein_clean.pdb"
        )

        if fallback.is_file():

            protein_path = fallback

        else:

            st.error(
                "The cleaned protein PDB "
                "could not be found."
            )

            st.code(
                str(protein_path)
            )

            return

    # --------------------------------------------------------
    # READ FILES
    # --------------------------------------------------------

    try:

        docking_text = (
            _read_text_file(
                output_path
            )
        )

        protein_text = (
            _read_text_file(
                protein_path
            )
        )

    except Exception as exc:

        st.error(
            "Could not read the "
            "docking files."
        )

        st.code(
            f"Docking output:\n"
            f"{output_path}\n\n"
            f"Protein PDB:\n"
            f"{protein_path}\n\n"
            f"Error:\n"
            f"{exc}"
        )

        return

    # --------------------------------------------------------
    # EXTRACT POSES
    # --------------------------------------------------------

    poses = (
        _extract_vina_poses(
            docking_text
        )
    )

    if not poses:

        st.error(
            "No Vina poses were found "
            "in the docking output."
        )

        return

    # --------------------------------------------------------
    # VIEWER CONTROLS
    # --------------------------------------------------------

    control_1, control_2, control_3 = (
        st.columns(3)
    )

    with control_1:

        pose_number = st.selectbox(
            "Docking pose",
            [
                pose["pose"]
                for pose in poses
            ],
            key="three_d_pose_selection",
        )

    with control_2:

        protein_style = st.selectbox(
            "Protein representation",
            [
                "Cartoon",
                "Surface",
                "Sticks",
            ],
            key="three_d_protein_style",
        )

    with control_3:

        ligand_style = st.selectbox(
            "Ligand representation",
            [
                "Sticks",
                "Ball & Stick",
                "Spheres",
            ],
            key="three_d_ligand_style",
        )

    dark_background = st.toggle(
        "Dark 3D viewer",
        value=True,
        key="three_d_dark_background",
    )

    selected_pose = next(
        pose
        for pose in poses
        if pose["pose"]
        == pose_number
    )

    # --------------------------------------------------------
    # CONVERT SELECTED VINA POSE
    # --------------------------------------------------------

    ligand_pdb = (
        pdbqt_pose_to_pdb(
            selected_pose["pdbqt"]
        )
    )

    if not ligand_pdb.strip():

        st.error(
            "The selected Vina pose "
            "contains no readable ligand atoms."
        )

        return

    # --------------------------------------------------------
    # METRICS
    # --------------------------------------------------------

    metric_1, metric_2, metric_3, metric_4 = (
        st.columns(4)
    )

    with metric_1:

        st.metric(
            "Protein",
            str(
                selected["Protein"]
            ),
        )

    with metric_2:

        st.metric(
            "Compound",
            str(
                selected["Compound"]
            ),
        )

    with metric_3:

        if (
            selected_pose["score"]
            is not None
        ):

            st.metric(
                "Binding affinity",
                (
                    f"{selected_pose['score']:.3f}"
                    " kcal/mol"
                ),
            )

        else:

            st.metric(
                "Binding affinity",
                "N/A",
            )

    with metric_4:

        st.metric(
            "Pose",
            str(
                pose_number
            ),
        )

    # --------------------------------------------------------
    # VIEWER + INFO
    # --------------------------------------------------------

    viewer_column, info_column = (
        st.columns(
            [3.4, 1.1],
            gap="large",
        )
    )

    with viewer_column:

        viewer = (
            _create_docking_viewer(
                protein_data=protein_text,
                ligand_data=ligand_pdb,
                protein_style=protein_style,
                ligand_style=ligand_style,
                dark_background=dark_background,
            )
        )

        components.html(
            viewer._make_html(),
            height=710,
            scrolling=False,
        )

    with info_column:

        st.markdown(
            "#### 🔬 Docking information"
        )

        st.write(
            f"**Protein:** "
            f"{selected['Protein']}"
        )

        st.write(
            f"**Compound:** "
            f"{selected['Compound']}"
        )

        st.write(
            f"**Chain:** "
            f"{selected.get('Chain', 'N/A')}"
        )

        st.write(
            f"**Pose:** "
            f"{pose_number}"
        )

        if (
            selected_pose["score"]
            is not None
        ):

            st.write(
                f"**Affinity:** "
                f"{selected_pose['score']:.3f} "
                "kcal/mol"
            )

        if (
            selected_pose["rmsd_lb"]
            is not None
        ):

            st.write(
                f"**RMSD LB:** "
                f"{selected_pose['rmsd_lb']:.3f} Å"
            )

        if (
            selected_pose["rmsd_ub"]
            is not None
        ):

            st.write(
                f"**RMSD UB:** "
                f"{selected_pose['rmsd_ub']:.3f} Å"
            )

        st.write(
            f"**PDB:** "
            f"{protein_path.name}"
        )

    # --------------------------------------------------------
    # POSE TABLE
    # --------------------------------------------------------

    pose_rows = []

    for pose in poses:

        pose_rows.append(
            {
                "Pose": pose["pose"],

                "Binding affinity (kcal/mol)": (
                    round(
                        pose["score"],
                        3,
                    )
                    if pose["score"]
                    is not None
                    else None
                ),

                "RMSD lower bound (Å)": (
                    round(
                        pose["rmsd_lb"],
                        3,
                    )
                    if pose["rmsd_lb"]
                    is not None
                    else None
                ),

                "RMSD upper bound (Å)": (
                    round(
                        pose["rmsd_ub"],
                        3,
                    )
                    if pose["rmsd_ub"]
                    is not None
                    else None
                ),
            }
        )

    st.markdown(
        "#### 📊 Vina poses"
    )

    pose_df = pd.DataFrame(
        pose_rows
    )

    st.dataframe(
        pose_df.style.format(
            {
                "Binding affinity (kcal/mol)": "{:.3f}",
                "RMSD lower bound (Å)": "{:.3f}",
                "RMSD upper bound (Å)": "{:.3f}",
            }
        ),
        use_container_width=True,
        hide_index=True,
    )

    # --------------------------------------------------------
    # FILE INFORMATION
    # --------------------------------------------------------

    with st.expander(
        "Docking files"
    ):

        st.write(
            f"**Protein:** "
            f"`{protein_path}`"
        )

        st.write(
            f"**Docking output:** "
            f"`{output_path}`"
        )


# ============================================================
# MAIN UI
# ============================================================

st.title(
    "🧬 DistDockNet"
)

st.caption(
    "Protein–compound molecular docking "
    "and redocking validation"
)


# ============================================================
# VINA CHECK
# ============================================================

if not VINA_EXE.exists():

    st.error(
        f"vina.exe was not found at:\n"
        f"{VINA_EXE}\n\n"
        "Place vina.exe in the same folder "
        "as app.py."
    )


# ============================================================
# SIDEBAR
# ============================================================

st.sidebar.header(
    "Docking settings"
)

workflow_options = [
    "Standard docking",
    "🧬 Disease-Guided Automated Docking",
    "Redocking validation",
    "🌐 Step 12: PPI & Network Analysis",
]

if st.session_state.get("jump_to_step12"):
    st.session_state["workflow_mode"] = "🌐 Step 12: PPI & Network Analysis"
    st.session_state["jump_to_step12"] = False

current_mode = st.session_state.get("workflow_mode", "Standard docking")
default_idx = (
    workflow_options.index(current_mode)
    if current_mode in workflow_options
    else 0
)

mode = st.sidebar.selectbox(
    "Workflow",
    workflow_options,
    index=default_idx,
)
st.session_state["workflow_mode"] = mode


site_mode = "Blind docking"


if mode == "Standard docking":

    site_mode = st.sidebar.selectbox(
        "Docking mode",
        [
            "Blind docking",
            "Reference ligand",
        ],
    )


chain_id = st.sidebar.text_input(
    "Preferred protein chain",
    value="A",
)


exhaustiveness = st.sidebar.slider(
    "Exhaustiveness",
    min_value=1,
    max_value=32,
    value=16,
)


num_modes = st.sidebar.slider(
    "Number of poses",
    min_value=1,
    max_value=20,
    value=10,
)


# ============================================================
# STEP 12: TARGET-COMPOUND INTERSECTION, PPI & FUNCTIONAL ANALYSIS
# ============================================================

if mode == "🌐 Step 12: PPI & Network Analysis":

    st.markdown(
        """
        <style>
        .step12-banner {
            background: linear-gradient(135deg, #0f172a 0%, #1e1b4b 50%, #0f172a 100%);
            border-left: 5px solid #38bdf8;
            border-radius: 12px;
            padding: 22px 28px;
            margin-bottom: 25px;
            box-shadow: 0 8px 24px rgba(0, 0, 0, 0.4);
        }
        .step12-title {
            font-size: 26px;
            font-weight: 700;
            color: #f8fafc;
            margin: 0 0 6px 0;
            letter-spacing: -0.5px;
        }
        .step12-subtitle {
            color: #94a3b8;
            font-size: 14px;
            margin: 0;
        }
        .step12-badge {
            background: rgba(56, 189, 248, 0.15);
            color: #38bdf8;
            border: 1px solid rgba(56, 189, 248, 0.35);
            padding: 3px 10px;
            border-radius: 20px;
            font-size: 11px;
            font-weight: 600;
            display: inline-block;
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        </style>
        <div class="step12-banner">
            <span class="step12-badge">Pipeline Architecture • Step 12</span>
            <div class="step12-title">🌐 Target–Compound Intersection, PPI & Functional Analysis</div>
            <p class="step12-subtitle">
                Systems pharmacology engine integrating STRING interactomes, master hub-target detection,
                KEGG & Gene Ontology functional enrichment, and bipartite polypharmacology scoring.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    if not _ANALYSIS_OK:
        st.error(
            "Analysis package failed to load. Please ensure networkx and pandas are installed in your environment."
        )
        st.stop()

    # ------------------------------------------------------------
    # Data Source & Configuration
    # ------------------------------------------------------------
    st.subheader("1. Input Data Source")

    default_batch_path = BASE_DIR / "results" / "batch_results.csv"
    has_batch = default_batch_path.exists()

    source_options = []
    if has_batch:
        source_options.append("Batch Docking Matrix (results/batch_results.csv)")
    source_options.extend([
        "SQLite Docking History (distdocknet.db)",
        "Upload Custom Docking Results CSV",
        "Target Genes Only (PPI & Enrichment without docking)",
    ])

    data_source = st.radio(
        "Select Data Source",
        source_options,
        index=0,
        horizontal=True,
    )

    docking_df = pd.DataFrame()
    target_genes = []

    c_cfg1, c_cfg2 = st.columns(2)

    with c_cfg1:
        string_score = st.slider(
            "STRING PPI Minimum Confidence Score",
            min_value=150,
            max_value=900,
            value=400,
            step=50,
            help="400 = Medium confidence, 700 = High confidence, 900 = Highest confidence.",
        )

    with c_cfg2:
        affinity_threshold = st.slider(
            "Strong Binding Affinity Cutoff (kcal/mol)",
            min_value=-14.0,
            max_value=-4.0,
            value=-7.0,
            step=0.5,
            help="Interactions at or below this energy threshold are treated as significant binding hits.",
        )

    if data_source == "Batch Docking Matrix (results/batch_results.csv)":
        try:
            docking_df = pd.read_csv(default_batch_path)
            st.success(f"✓ Loaded {len(docking_df)} batch docking records from `{default_batch_path.name}`.")
            if "protein" in docking_df.columns:
                target_genes = [
                    normalize_target_to_gene(p)
                    for p in docking_df["protein"].dropna().unique()
                ]
            with st.expander("Preview Batch Docking Table", expanded=False):
                st.dataframe(docking_df.head(15), use_container_width=True)
        except Exception as e:
            st.error(f"Error reading batch results: {e}")

    elif data_source == "SQLite Docking History (distdocknet.db)":
        from database import get_docking_history
        rows, cols = get_docking_history(limit=500)
        if rows:
            docking_df = pd.DataFrame(rows, columns=cols)
            # Filter successful runs
            if "status" in docking_df.columns:
                docking_df = docking_df[docking_df["status"] == "SUCCESS"]
            st.success(f"✓ Loaded {len(docking_df)} successful docking runs from SQLite database.")
            if "protein_name" in docking_df.columns:
                target_genes = [
                    normalize_target_to_gene(p)
                    for p in docking_df["protein_name"].dropna().unique()
                ]
            with st.expander("Preview Database History", expanded=False):
                st.dataframe(docking_df.head(15), use_container_width=True)
        else:
            st.warning("No docking runs found in database. Run docking first or select another data source.")

    elif data_source == "Upload Custom Docking Results CSV":
        uploaded_csv = st.file_uploader(
            "Upload Docking CSV (must have columns for protein, compound, and affinity)",
            type=["csv"],
        )
        if uploaded_csv is not None:
            try:
                docking_df = pd.read_csv(uploaded_csv)
                st.success(f"✓ Uploaded {len(docking_df)} records.")
                # Auto-detect protein column
                for col in ["protein", "protein_name", "target", "gene"]:
                    if col in docking_df.columns:
                        target_genes = [
                            normalize_target_to_gene(p)
                            for p in docking_df[col].dropna().unique()
                        ]
                        break
                with st.expander("Preview Uploaded CSV", expanded=False):
                    st.dataframe(docking_df.head(15), use_container_width=True)
            except Exception as e:
                st.error(f"Error reading CSV: {e}")

    else:
        # Target Genes Only
        gene_input = st.text_area(
            "Enter Disease Target Genes (comma-separated)",
            value="INS, AKT1, PPARG, TNF, IL6, MTOR, MAPK1, TP53, VEGFA, EGFR, MMP13, STAT3",
            help="Standard HGNC gene symbols for STRING interactome & functional enrichment.",
        )
        target_genes = [
            normalize_target_to_gene(g.strip())
            for g in gene_input.split(",")
            if g.strip()
        ]

    # De-duplicate target genes
    target_genes = sorted(list(set(target_genes)))

    st.markdown("---")

    # ------------------------------------------------------------
    # Run Analysis Button
    # ------------------------------------------------------------
    run_analysis = st.button("🚀 Run Step 12 Systems Analysis", type="primary", use_container_width=True)

    if run_analysis or "step12_results" in st.session_state:

        if run_analysis:
            if not target_genes:
                st.warning("No target genes identified to analyze. Please provide target genes or docking records.")
                st.stop()

            with st.spinner("Querying STRING DB & constructing Systems Biology interactome..."):
                # 1. Fetch PPI and Centralities
                ppi_edges = fetch_string_ppi(
                    genes=target_genes,
                    species=9606,
                    required_score=string_score,
                )
                G_ppi = build_ppi_graph(ppi_edges, user_genes=target_genes)
                centralities_df = calculate_network_centralities(G_ppi)

                # 2. Fetch Functional Enrichment (KEGG, GO, Reactome)
                enrichment_dfs = fetch_functional_enrichment(
                    genes=target_genes,
                    species=9606,
                    fdr_threshold=0.05,
                )

                # 3. Target-Compound Intersection
                G_tc = None
                poly_df = pd.DataFrame()
                if not docking_df.empty:
                    # Detect column names
                    t_col = "protein" if "protein" in docking_df.columns else "protein_name"
                    c_col = "compound" if "compound" in docking_df.columns else "compound_name"
                    a_col = "binding_affinity"

                    G_tc = build_target_compound_network(
                        docking_df=docking_df,
                        affinity_threshold=affinity_threshold,
                        target_col=t_col,
                        compound_col=c_col,
                        affinity_col=a_col,
                    )
                    poly_df = rank_compounds_by_polypharmacology(
                        docking_df=docking_df,
                        hub_df=centralities_df,
                        affinity_threshold=affinity_threshold,
                        target_col=t_col,
                        compound_col=c_col,
                        affinity_col=a_col,
                    )

                st.session_state["step12_results"] = {
                    "target_genes": target_genes,
                    "ppi_edges": ppi_edges,
                    "G_ppi": G_ppi,
                    "centralities_df": centralities_df,
                    "enrichment_dfs": enrichment_dfs,
                    "G_tc": G_tc,
                    "poly_df": poly_df,
                    "has_docking": not docking_df.empty,
                }

        # Retrieve stored results
        results = st.session_state.get("step12_results", {})
        G_ppi = results.get("G_ppi")
        centralities_df = results.get("centralities_df", pd.DataFrame())
        enrichment_dfs = results.get("enrichment_dfs", {})
        G_tc = results.get("G_tc")
        poly_df = results.get("poly_df", pd.DataFrame())
        has_docking = results.get("has_docking", False)

        st.subheader("2. Analysis Results & Interactome Dashboards")

        tab_ppi, tab_enrich, tab_poly = st.tabs([
            "🕸️ Protein–Protein Interaction (PPI) & Hub Targets",
            "🧬 Functional & Pathway Enrichment (KEGG & GO)",
            "💊 Target–Compound Intersection & Polypharmacology",
        ])

        # ============================================================
        # TAB 1: PPI NETWORK & MASTER HUBS
        # ============================================================
        with tab_ppi:
            st.markdown("### 🕸️ Protein–Protein Interaction Interactome")
            st.caption(
                "Direct and functional protein interactions derived from STRING DB. Master hub targets are detected via graph topology and centrality metrics."
            )

            # Metrics row
            m1, m2, m3, m4 = st.columns(4)
            num_nodes = len(G_ppi.nodes) if G_ppi else 0
            num_edges = len(G_ppi.edges) if G_ppi else 0
            hub_count = len(centralities_df[centralities_df["is_hub"]]) if not centralities_df.empty else 0
            avg_deg = (2 * num_edges / num_nodes) if num_nodes > 0 else 0

            with m1:
                st.metric("Target Genes", num_nodes)
            with m2:
                st.metric("Interactome Edges", num_edges)
            with m3:
                st.metric("Average Degree", f"{avg_deg:.1f}")
            with m4:
                st.metric("Master Hub Targets", hub_count)

            # Vis.js interactive network viewer
            if G_ppi and num_nodes > 0:
                html_ppi = generate_network_html(
                    G_ppi,
                    title="Target PPI Interactome (STRING)",
                    network_type="ppi",
                    height_px=560,
                )
                components.html(html_ppi, height=580, scrolling=False)
            else:
                st.info("No PPI interactions found at this confidence cutoff. Try lowering the STRING confidence score.")

            st.markdown("#### 🏆 Master Hub Targets & Topological Centralities")
            if not centralities_df.empty:
                st.dataframe(
                    centralities_df.style.highlight_max(
                        subset=["degree", "hub_score", "betweenness_centrality"],
                        color="#1e3a8a",
                    ),
                    use_container_width=True,
                    hide_index=True,
                )
                st.download_button(
                    "📥 Download Centralities & Hub Targets (CSV)",
                    data=centralities_df.to_csv(index=False).encode("utf-8"),
                    file_name="distdocknet_hub_targets_centralities.csv",
                    mime="text/csv",
                )

        # ============================================================
        # TAB 2: FUNCTIONAL & PATHWAY ENRICHMENT
        # ============================================================
        with tab_enrich:
            st.markdown("### 🧬 Functional & Pathway Enrichment Analysis")
            st.caption(
                "Identifies biologically overrepresented pathways and Gene Ontology terms among target proteins."
            )

            available_categories = [
                k for k, v in enrichment_dfs.items() if not v.empty
            ]

            if available_categories:
                selected_cat = st.radio(
                    "Select Enrichment Domain",
                    available_categories,
                    horizontal=True,
                )

                cat_df = enrichment_dfs.get(selected_cat, pd.DataFrame())

                if not cat_df.empty:
                    chart = get_enrichment_altair_chart(
                        cat_df,
                        title=f"Top Enriched {selected_cat} (-log10 FDR)",
                        top_n=10,
                    )
                    if chart:
                        st.altair_chart(chart, use_container_width=True)

                    st.markdown(f"#### All Enriched {selected_cat} Terms ({len(cat_df)})")
                    st.dataframe(
                        cat_df,
                        use_container_width=True,
                        hide_index=True,
                    )
                    st.download_button(
                        f"📥 Download {selected_cat} Results (CSV)",
                        data=cat_df.to_csv(index=False).encode("utf-8"),
                        file_name=f"distdocknet_{selected_cat.lower().replace(' ', '_')}.csv",
                        mime="text/csv",
                    )
            else:
                st.info("No statistically significant enriched pathways found (FDR <= 0.05).")

        # ============================================================
        # TAB 3: TARGET–COMPOUND INTERSECTION & POLYPHARMACOLOGY
        # ============================================================
        with tab_poly:
            st.markdown("### 💊 Target–Compound Intersection & Polypharmacology")
            st.caption(
                "Bipartite interaction network mapping nutraceutical compounds against target proteins with polypharmacology scoring prioritizing compounds hitting multiple disease hub targets."
            )

            if has_docking and not poly_df.empty:
                p1, p2, p3, p4 = st.columns(4)
                multi_hits = len(poly_df[poly_df["strong_targets_count"] >= 2])
                top_lead = poly_df.iloc[0]["compound"] if len(poly_df) > 0 else "N/A"
                best_aff = poly_df["best_affinity_kcal_mol"].min() if len(poly_df) > 0 else 0.0

                with p1:
                    st.metric("Tested Compounds", len(poly_df))
                with p2:
                    st.metric("Multi-Target Leads (≥2)", multi_hits)
                with p3:
                    st.metric("Top Ranked Compound", top_lead)
                with p4:
                    st.metric("Best Affinity", f"{best_aff:.2f} kcal/mol")

                # Bipartite Network Visualizer
                if G_tc and len(G_tc.nodes) > 0:
                    html_tc = generate_network_html(
                        G_tc,
                        title=f"Target–Compound Interaction Network (Affinity ≤ {affinity_threshold} kcal/mol)",
                        network_type="target_compound",
                        height_px=560,
                    )
                    components.html(html_tc, height=580, scrolling=False)

                st.markdown("#### 🏆 Multi-Target Polypharmacology Ranking")
                st.dataframe(
                    poly_df.style.highlight_max(
                        subset=["poly_score", "strong_targets_count"],
                        color="#064e3b",
                    ),
                    use_container_width=True,
                    hide_index=True,
                )
                st.download_button(
                    "📥 Download Polypharmacology Rankings (CSV)",
                    data=poly_df.to_csv(index=False).encode("utf-8"),
                    file_name="distdocknet_polypharmacology_rankings.csv",
                    mime="text/csv",
                )
            else:
                st.info(
                    "Docking results are required for Target–Compound intersection analysis. Select 'Batch Docking Matrix' or 'Upload Custom Docking CSV' as your input data source above."
                )


# ============================================================
# DISEASE-GUIDED AUTOMATED DOCKING WORKFLOW
# ============================================================

elif mode == "🧬 Disease-Guided Automated Docking":

    from disease_docking_ui import render_disease_docking_workflow

    render_disease_docking_workflow(
        chain_id=chain_id,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        vina_exe=VINA_EXE,
        fetch_pdb_func=fetch_pdb,
        run_standard_docking_func=run_standard_docking,
        render_3d_viewer_func=_render_docking_3d_viewer,
    )


# ============================================================
# REDOCKING WORKFLOW
# ============================================================

elif mode == "Redocking validation":

    st.info(
        "Redocking validation extracts the known ligand "
        "from the experimental PDB, redocks it, and "
        "compares the predicted pose with the experimental pose."
    )

    st.subheader(
        "Redocking validation input"
    )

    c1, c2 = st.columns(
        2
    )

    with c1:

        protein_name = st.text_input(
            "Protein name",
            value="MMP-13",
        )

        pdb_id = st.text_input(
            "PDB ID",
            value="2OW9",
        )

    with c2:

        reference_ligand = st.text_input(
            "Reference ligand CCD code",
            value="SP6",
        )

        st.write(
            "Expected example: "
            "MMP-13 / 2OW9 / SP6"
        )

    st.divider()

    if st.button(
        "🧪 Run redocking validation",
        type="primary",
        use_container_width=True,
    ):

        if not VINA_EXE.exists():

            st.error(
                "vina.exe is missing."
            )

            st.stop()

        if (
            not protein_name.strip()
            or not pdb_id.strip()
            or not reference_ligand.strip()
        ):

            st.error(
                "Enter protein name, "
                "PDB ID, and reference "
                "ligand CCD code."
            )

            st.stop()

        try:

            with st.spinner(
                f"Getting "
                f"{pdb_id.upper()} "
                f"from RCSB..."
            ):

                pdb_file = fetch_pdb(
                    pdb_id
                )

            with st.spinner(
                f"Redocking "
                f"{reference_ligand.upper()} "
                f"into "
                f"{protein_name}..."
            ):

                result = (
                    run_redocking_validation(
                        protein_name=protein_name,
                        pdb_file=pdb_file,
                        reference_ligand=reference_ligand,
                        chain_id=(
                            chain_id.strip()
                            or "A"
                        ),
                        exhaustiveness=exhaustiveness,
                        num_modes=num_modes,
                    )
                )

            st.success(
                "Redocking completed successfully."
            )

            c1, c2, c3 = st.columns(
                3
            )

            with c1:

                st.metric(
                    "Vina score",
                    (
                        f"{result['binding_affinity']:.3f}"
                        " kcal/mol"
                    ),
                )

            with c2:

                st.metric(
                    "Pose RMSD",
                    (
                        f"{result['rmsd']:.3f}"
                        " Å"
                    ),
                )

            with c3:

                st.metric(
                    "Matched atoms",
                    str(
                        result[
                            "matched_atoms"
                        ]
                    ),
                )

            if (
                result["rmsd"]
                <= 2.0
            ):

                st.success(
                    "RMSD ≤ 2.0 Å: "
                    "the redocked pose reproduces "
                    "the experimental ligand position "
                    "reasonably well."
                )

            else:

                st.warning(
                    "RMSD > 2.0 Å: "
                    "the redocked pose differs "
                    "substantially from the experimental pose."
                )

            st.write(
                "**RMSD calculation:**",
                result[
                    "rmsd_method"
                ],
            )

            st.write(
                "**Experimental ligand:**",
                result[
                    "experimental_ligand"
                ],
            )

            st.write(
                "**Predicted docking output:**",
                result[
                    "predicted_output"
                ],
            )

            st.write(
                "**Docking box:**",
                result[
                    "box"
                ],
            )

            st.session_state[
                "redocking_result"
            ] = result

        except Exception as exc:

            st.error(
                "Redocking validation failed."
            )

            st.exception(
                exc
            )

    # --------------------------------------------------------
    # REDOCKING FILES
    # --------------------------------------------------------

    if (
        "redocking_result"
        in st.session_state
    ):

        result = (
            st.session_state[
                "redocking_result"
            ]
        )

        st.divider()

        st.subheader(
            "Redocking files"
        )

        file_rows = pd.DataFrame(
            [
                {
                    "File":
                        "Experimental ligand",
                    "Path":
                        result[
                            "experimental_ligand"
                        ],
                },
                {
                    "File":
                        "Prepared reference ligand",
                    "Path":
                        result[
                            "reference_ligand_pdbqt"
                        ],
                },
                {
                    "File":
                        "Clean protein",
                    "Path":
                        result[
                            "protein_pdb"
                        ],
                },
                {
                    "File":
                        "Receptor PDBQT",
                    "Path":
                        result[
                            "receptor_pdbqt"
                        ],
                },
                {
                    "File":
                        "Predicted pose",
                    "Path":
                        result[
                            "predicted_output"
                        ],
                },
            ]
        )

        st.dataframe(
            file_rows,
            use_container_width=True,
            hide_index=True,
        )


# ============================================================
# STANDARD DOCKING WORKFLOW
# ============================================================

else:

    st.subheader(
        "1. Proteins"
    )

    protein_count = st.number_input(
        "Number of proteins",
        min_value=1,
        max_value=20,
        value=1,
        step=1,
    )

    proteins = []

    for i in range(
        int(protein_count)
    ):

        st.markdown(
            f"**Protein {i + 1}**"
        )

        c1, c2 = st.columns(
            2
        )

        with c1:

            name = st.text_input(
                "Protein name",
                value=(
                    f"Protein_{i + 1}"
                ),
                key=(
                    f"protein_name_{i}"
                ),
            )

        with c2:

            source = st.radio(
                "Protein source",
                [
                    "RCSB PDB",
                    "Upload PDB file",
                ],
                index=0,
                horizontal=True,
                key=(
                    f"protein_source_{i}"
                ),
            )

        pdb_id = ""

        uploaded = None

        if (
            source
            == "RCSB PDB"
        ):

            pdb_id = st.text_input(
                "PDB ID",
                placeholder="Example: 2OW9",
                key=(
                    f"pdb_id_{i}"
                ),
            )

        else:

            uploaded = st.file_uploader(
                "Upload PDB file",
                type=[
                    "pdb",
                    "ent",
                ],
                key=(
                    f"protein_file_{i}"
                ),
                help=(
                    "Upload a local "
                    ".pdb or .ent structure."
                ),
            )

        proteins.append(
            {
                "name": name,
                "pdb_id": pdb_id,
                "uploaded": uploaded,
                "source": source,
            }
        )

    # ========================================================
    # COMPOUNDS
    # ========================================================

    st.subheader(
        "2. Compounds"
    )

    compound_count = st.number_input(
        "Number of compounds",
        min_value=1,
        max_value=20,
        value=1,
        step=1,
    )

    compounds = []

    for i in range(
        int(compound_count)
    ):

        st.markdown(
            f"**Compound {i + 1}**"
        )

        name = st.text_input(
            "Compound name",
            value=(
                f"Compound_{i + 1}"
            ),
            key=(
                f"compound_name_{i}"
            ),
        )

        smiles = st.text_area(
            "SMILES",
            placeholder=(
                "Paste the compound "
                "SMILES here"
            ),
            key=(
                f"compound_smiles_{i}"
            ),
            height=70,
        )

        compounds.append(
            {
                "name": name,
                "smiles": smiles,
            }
        )

    # ========================================================
    # REFERENCE LIGAND MODE
    # ========================================================

    reference_ligand = ""

    if (
        site_mode
        == "Reference ligand"
    ):

        reference_ligand = st.text_input(
            "Reference ligand CCD code",
            placeholder="Example: SP6",
        )

    # ========================================================
    # RUN DOCKING
    # ========================================================

    st.divider()

    st.subheader(
        "3. Run docking"
    )

    if st.button(
        "🚀 Start docking",
        type="primary",
        use_container_width=True,
    ):

        if not VINA_EXE.exists():

            st.error(
                "vina.exe is missing."
            )

            st.stop()

        valid_compounds = [
            compound
            for compound in compounds
            if (
                compound["name"].strip()
                and compound["smiles"].strip()
            )
        ]

        if not valid_compounds:

            st.error(
                "Enter at least one "
                "valid compound."
            )

            st.stop()

        if (
            site_mode
            == "Reference ligand"
            and not reference_ligand.strip()
        ):

            st.error(
                "Enter a reference "
                "ligand CCD code."
            )

            st.stop()

        resolved_proteins = []

        try:

            for protein in proteins:

                if not protein[
                    "name"
                ].strip():

                    raise ValueError(
                        "Every protein "
                        "needs a name."
                    )

                # ------------------------------------------------
                # UPLOAD PDB
                # ------------------------------------------------

                if (
                    protein["source"]
                    == "Upload PDB file"
                ):

                    if (
                        protein["uploaded"]
                        is None
                    ):

                        raise ValueError(
                            f"Upload a PDB "
                            f"file for "
                            f"'{protein['name']}'."
                        )

                    destination = (
                        PROTEIN_CACHE
                        / safe_name(
                            protein["name"]
                        )
                        / "input.pdb"
                    )

                    destination.parent.mkdir(
                        parents=True,
                        exist_ok=True,
                    )

                    destination.write_bytes(
                        protein[
                            "uploaded"
                        ].getvalue()
                    )

                    pdb_file = destination

                # ------------------------------------------------
                # RCSB PDB
                # ------------------------------------------------

                else:

                    if not protein[
                        "pdb_id"
                    ].strip():

                        raise ValueError(
                            f"Enter a PDB ID "
                            f"for "
                            f"'{protein['name']}'."
                        )

                    with st.spinner(
                        f"Getting "
                        f"{protein['pdb_id'].upper()}"
                        f"..."
                    ):

                        pdb_file = fetch_pdb(
                            protein[
                                "pdb_id"
                            ]
                        )

                resolved_proteins.append(
                    {
                        "name":
                            protein[
                                "name"
                            ],
                        "pdb_file":
                            pdb_file,
                    }
                )

        except Exception as exc:

            st.error(
                str(exc)
            )

            st.stop()

        # --------------------------------------------------------
        # CREATE ALL JOBS
        # --------------------------------------------------------

        jobs = list(
            product(
                resolved_proteins,
                valid_compounds,
            )
        )

        st.info(
            f"Running "
            f"{len(jobs)} "
            f"docking combination(s): "
            f"{len(resolved_proteins)} "
            f"protein(s) × "
            f"{len(valid_compounds)} "
            f"compound(s)."
        )

        progress = st.progress(
            0
        )

        status_box = st.empty()

        results = []

        # ========================================================
        # DOCK EACH JOB
        # ========================================================

        for index, (
            protein,
            compound,
        ) in enumerate(
            jobs,
            start=1,
        ):

            status_box.write(
                f"Docking "
                f"{index}/{len(jobs)}: "
                f"{protein['name']} + "
                f"{compound['name']}"
            )

            try:

                result = (
                    run_standard_docking(
                        protein_name=
                            protein[
                                "name"
                            ],
                        pdb_file=
                            protein[
                                "pdb_file"
                            ],
                        compound_name=
                            compound[
                                "name"
                            ],
                        smiles=
                            compound[
                                "smiles"
                            ],
                        chain_id=(
                            chain_id.strip()
                            or "A"
                        ),
                        site_mode=(
                            "reference"
                            if (
                                site_mode
                                == "Reference ligand"
                            )
                            else "blind"
                        ),
                        reference_ligand=
                            reference_ligand,
                        exhaustiveness=
                            exhaustiveness,
                        num_modes=
                            num_modes,
                    )
                )

                results.append(
                    {
                        "Protein":
                            result[
                                "protein"
                            ],

                        "Compound":
                            result[
                                "compound"
                            ],

                        "Binding affinity (kcal/mol)":
                            result[
                                "binding_affinity"
                            ],

                        "Chain":
                            result[
                                "selected_chain"
                            ],

                        "Docking output":
                            result[
                                "output"
                            ],

                        "Protein PDB":
                            result[
                                "protein_pdb"
                            ],

                        "Ligand PDBQT":
                            result[
                                "ligand_pdbqt"
                            ],
                    }
                )

            except Exception as exc:

                results.append(
                    {
                        "Protein":
                            protein[
                                "name"
                            ],

                        "Compound":
                            compound[
                                "name"
                            ],

                        "Binding affinity (kcal/mol)":
                            None,

                        "Chain":
                            None,

                        "Docking output":
                            "",

                        "Protein PDB":
                            "",

                        "Ligand PDBQT":
                            "",

                        "Error":
                            str(exc),
                    }
                )

            progress.progress(
                index
                / len(jobs)
            )

        status_box.success(
            "Docking batch completed."
        )

        # ========================================================
        # RESULTS DATAFRAME
        # ========================================================

        results_df = pd.DataFrame(
            results
        )

        output_csv = (
            RESULTS_DIR
            / "ui_results.csv"
        )

        results_df.to_csv(
            output_csv,
            index=False,
        )

        st.session_state[
            "results_df"
        ] = results_df

        st.success(
            f"Completed "
            f"{len(jobs)} "
            f"docking combination(s)."
        )

    # ============================================================
    # DISPLAY RESULTS
    # ============================================================

    if (
        "results_df"
        in st.session_state
    ):

        st.divider()

        st.subheader(
            "4. Docking results"
        )

        results_df = (
            st.session_state[
                "results_df"
            ]
        )

        affinity = pd.to_numeric(
            results_df[
                "Binding affinity (kcal/mol)"
            ],
            errors="coerce",
        )

        successful = (
            results_df.loc[
                affinity.notna()
            ].copy()
        )

        # --------------------------------------------------------
        # SUCCESSFUL RESULTS
        # --------------------------------------------------------

        if not successful.empty:

            successful = (
                successful.sort_values(
                    "Binding affinity (kcal/mol)"
                )
            )

            # Explicitly display 3 decimal places.
            st.dataframe(
                successful.style.format(
                    {
                        "Binding affinity (kcal/mol)":
                            "{:.3f}",
                    }
                ),
                use_container_width=True,
                hide_index=True,
            )

            st.subheader(
                "Best docking score"
            )

            best = (
                successful.iloc[0]
            )

            best_score = float(
                best[
                    "Binding affinity (kcal/mol)"
                ]
            )

            st.metric(
                label=(
                    f"{best['Protein']} + "
                    f"{best['Compound']}"
                ),
                value=(
                    f"{best_score:.3f} "
                    "kcal/mol"
                ),
            )

        # --------------------------------------------------------
        # FAILED RESULTS
        # --------------------------------------------------------

        failed = (
            results_df.loc[
                affinity.isna()
            ]
        )

        if not failed.empty:

            st.warning(
                f"{len(failed)} "
                f"docking job(s) failed."
            )

            st.dataframe(
                failed,
                use_container_width=True,
                hide_index=True,
            )

        # --------------------------------------------------------
        # DOWNLOAD CSV
        # --------------------------------------------------------

        csv_data = (
            results_df
            .to_csv(
                index=False
            )
            .encode(
                "utf-8"
            )
        )

        st.download_button(
            "Download results CSV",
            data=csv_data,
            file_name=(
                "distdocknet_ui_results.csv"
            ),
            mime="text/csv",
        )

        # ========================================================
        # INTERACTIVE 3D VIEWER
        # ========================================================

        _render_docking_3d_viewer(
            results_df
        )


# ============================================================
# FOOTER
# ============================================================

st.divider()

st.caption(
    "Scores are Vina docking predictions and "
    "should not be interpreted as experimental "
    "binding affinities or pharmacological potency."
)