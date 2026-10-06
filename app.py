import io

import re

import subprocess

from pathlib import Path

from itertools import product



import numpy as np

import pandas as pd

import streamlit as st

from rdkit import Chem

from rdkit.Chem import AllChem

from Bio.PDB import PDBParser, PDBIO, Select

from meeko import MoleculePreparation, PDBQTWriterLegacy
from database import initialize_database, create_working_pdb_file


def smiles_to_sdf(smiles, output_file):
    print("\n[1] Converting SMILES to 3D structure...")

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Invalid SMILES string.")

    print("SMILES parsed successfully.")

    mol = Chem.AddHs(mol)

    print("Generating 3D coordinates...")
    result = AllChem.EmbedMolecule(mol, randomSeed=42)

    if result != 0:
        raise RuntimeError("Could not generate 3D coordinates for the compound.")

    print("3D coordinates generated.")

    print("Optimizing molecular geometry...")
    result = AllChem.UFFOptimizeMolecule(mol)

    if result != 0:
        print("Warning: UFF optimization did not fully converge.")

    writer = Chem.SDWriter(str(output_file))
    writer.write(mol)
    writer.close()

    print("3D SDF saved:", output_file)

    return output_file
BASE_DIR = Path(__file__).resolve().parent
VINA_EXE = BASE_DIR / "vina.exe"
RESULTS_DIR = BASE_DIR / "results" / "ui_runs"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# SQLite is the permanent protein store.
initialize_database()



STANDARD_AMINO_ACIDS = {

    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",

    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",

    "THR", "TRP", "TYR", "VAL"

}





class ProteinSelect(Select):

    def __init__(self, chain_id):

        self.chain_id = chain_id



    def accept_chain(self, chain):

        return chain.id == self.chain_id



    def accept_residue(self, residue):

        if residue.resname.strip() not in STANDARD_AMINO_ACIDS:

            return 0

        if residue.id[1] < 0:

            return 0

        return 1





def safe_name(value):

    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())

    return value[:80] or "item"





def fetch_pdb(pdb_id, working_directory):
    """Get a protein from SQLite/RCSB and create a working PDB file."""

    pdb_id = pdb_id.strip().upper()

    if not re.fullmatch(r"[0-9A-Za-z]{4}", pdb_id):
        raise ValueError("PDB ID must contain exactly 4 letters/numbers.")

    return create_working_pdb_file(pdb_id, working_directory)


def extract_and_clean_protein(input_pdb, output_pdb, chain_id="A"):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", str(input_pdb))

    model = next(structure.get_models(), None)

    if model is None:

        raise ValueError("Protein has no model.")



    available = {}

    for chain in model:

        count = sum(

            1 for residue in chain

            if residue.resname.strip() in STANDARD_AMINO_ACIDS

        )

        if count:

            available[chain.id] = count



    if not available:

        raise ValueError("No protein chains found.")



    selected = chain_id if chain_id in available else max(available, key=available.get)



    io_obj = PDBIO()

    io_obj.set_structure(structure)

    io_obj.save(str(output_pdb), ProteinSelect(selected))

    return selected, available





def extract_reference_ligand(input_pdb, output_pdb, ligand_name, chain_id=None):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("complex", str(input_pdb))



    matches = []

    for model in structure:

        for chain in model:

            if chain_id and chain.id != chain_id:

                continue

            for residue in chain:

                if residue.resname.strip() == ligand_name.strip():

                    matches.append((chain.id, residue))



    if not matches:

        raise ValueError(

            f"Reference ligand '{ligand_name}' was not found in the PDB."

        )



    selected_chain, selected_residue = matches[0]



    class SingleLigandSelect(Select):

        def accept_model(self, model):

            return model.id == 0



        def accept_chain(self, chain):

            return chain.id == selected_chain



        def accept_residue(self, residue):

            return (

                residue.resname.strip() == ligand_name.strip()

                and residue.id == selected_residue.id

            )



    io_obj = PDBIO()

    io_obj.set_structure(structure)

    io_obj.save(str(output_pdb), SingleLigandSelect())



    return selected_chain





def load_reference_ligand_as_rdkit(ligand_pdb):

    mol = Chem.MolFromPDBFile(

        str(ligand_pdb),

        removeHs=False,

        sanitize=False,

    )

    if mol is None:

        raise ValueError("RDKit could not read the extracted reference ligand.")



    try:

        Chem.SanitizeMol(mol)

    except Exception:

        pass



    if mol.GetNumConformers() == 0:

        raise ValueError("Reference ligand has no 3D coordinates.")



    return mol





def prepare_rdkit_ligand(mol, output_pdbqt):

    """

    Prepare an RDKit ligand molecule for AutoDock Vina using Meeko.



    Meeko requires explicit hydrogens, so we explicitly add H atoms

    before ligand preparation.

    """

    if mol is None:

        raise ValueError("Cannot prepare ligand because the RDKit molecule is None.")



    print("\nPreparing ligand with Meeko...")



    # Make a writable copy so the original molecule is not modified.

    mol = Chem.Mol(mol)



    # Meeko requires explicit hydrogens.

    try:

        mol = Chem.AddHs(mol)

    except Exception as e:

        raise RuntimeError(

            f"Could not add explicit hydrogens to the ligand: {e}"

        ) from e



    print("Explicit hydrogens added.")



    # Check that the molecule now contains explicit H atoms.

    hydrogen_count = sum(

        1 for atom in mol.GetAtoms()

        if atom.GetAtomicNum() == 1

    )



    print(f"Explicit hydrogen atoms: {hydrogen_count}")



    if hydrogen_count == 0:

        raise RuntimeError(

            "No explicit hydrogen atoms were added. "

            "Meeko requires explicit hydrogens."

        )



    # Prepare molecule with Meeko.

    print("Running Meeko ligand preparation...")



    preparator = MoleculePreparation()



    try:

        setups = preparator.prepare(mol)

    except Exception as e:

        raise RuntimeError(

            f"Meeko ligand preparation failed: {e}"

        ) from e



    if not setups:

        raise RuntimeError(

            "Meeko did not generate a ligand setup."

        )



    # Convert Meeko setup to PDBQT.

    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(

        setups[0]

    )



    if not is_ok:

        raise RuntimeError(

            f"Meeko PDBQT writing failed: {error_msg}"

        )



    output_pdbqt = str(output_pdbqt)



    with open(output_pdbqt, "w") as f:

        f.write(pdbqt_string)



    print("Ligand preparation completed successfully.")

    print("Saved as:", output_pdbqt)



    return output_pdbqt





def prepare_ligand(sdf_file, output_pdbqt):

    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)

    mol = supplier[0]

    if mol is None:

        raise ValueError("Could not read generated ligand structure.")

    return prepare_rdkit_ligand(mol, output_pdbqt)





def prepare_receptor(pdb_file, output_prefix):

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



    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:

        details = result.stderr or result.stdout or "Unknown Meeko error."

        raise RuntimeError(

            "Meeko receptor preparation failed.\n\n" + details

        )



    receptor = Path(str(output_prefix) + ".pdbqt")

    if not receptor.exists():

        raise RuntimeError("Meeko completed but receptor PDBQT was not created.")

    return receptor





def protein_coordinates(pdb_file):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", str(pdb_file))

    coordinates = []



    for model in structure:

        for chain in model:

            for residue in chain:

                if residue.resname.strip() not in STANDARD_AMINO_ACIDS:

                    continue

                for atom in residue:

                    coordinates.append(atom.coord)



    if not coordinates:

        raise ValueError("No protein coordinates found.")

    return coordinates





def calculate_blind_box(pdb_file):

    coordinates = protein_coordinates(pdb_file)

    xs = [float(c[0]) for c in coordinates]

    ys = [float(c[1]) for c in coordinates]

    zs = [float(c[2]) for c in coordinates]



    padding = 5.0

    return {

        "center_x": (min(xs) + max(xs)) / 2,

        "center_y": (min(ys) + max(ys)) / 2,

        "center_z": (min(zs) + max(zs)) / 2,

        "size_x": (max(xs) - min(xs)) + 2 * padding,

        "size_y": (max(ys) - min(ys)) + 2 * padding,

        "size_z": (max(zs) - min(zs)) + 2 * padding,

    }





def calculate_reference_ligand_box(pdb_file, ligand_resname, padding=8.0, chain_id=None):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("protein", str(pdb_file))

    coordinates = []



    for model in structure:

        for chain in model:

            if chain_id and chain.id != chain_id:

                continue

            for residue in chain:

                if residue.resname.strip() == ligand_resname.strip():

                    for atom in residue:

                        coordinates.append(atom.coord)



    if not coordinates:

        raise ValueError(

            f"Reference ligand '{ligand_resname}' was not found."

        )



    xs = [float(c[0]) for c in coordinates]

    ys = [float(c[1]) for c in coordinates]

    zs = [float(c[2]) for c in coordinates]



    return {

        "center_x": (min(xs) + max(xs)) / 2,

        "center_y": (min(ys) + max(ys)) / 2,

        "center_z": (min(zs) + max(zs)) / 2,

        "size_x": (max(xs) - min(xs)) + 2 * padding,

        "size_y": (max(ys) - min(ys)) + 2 * padding,

        "size_z": (max(zs) - min(zs)) + 2 * padding,

    }





def run_vina(receptor, ligand, box, output_file, exhaustiveness=16, num_modes=10):

    config_file = output_file.parent / "vina_config.txt"

    config_text = f"""receptor = {receptor}

ligand = {ligand}



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

    config_file.write_text(config_text, encoding="utf-8")



    command = [

        str(VINA_EXE),

        "--config",

        str(config_file),

        "--out",

        str(output_file),

    ]



    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:

        details = result.stderr or result.stdout or "Unknown Vina error."

        raise RuntimeError("Vina docking failed.\n\n" + details)



    return result.stdout





def parse_vina_score(pdbqt_file):

    pattern = re.compile(r"REMARK VINA RESULT:\s+(-?\d+(?:\\.\d+)?)")

    for line in Path(pdbqt_file).read_text(

        encoding="utf-8", errors="ignore"

    ).splitlines():

        match = pattern.search(line)

        if match:

            return float(match.group(1))

    raise ValueError("Could not find VINA RESULT in docking output.")





def parse_pdbqt_atoms(pdbqt_file):

    atoms = []

    for line in Path(pdbqt_file).read_text(

        encoding="utf-8", errors="ignore"

    ).splitlines():

        if line.startswith(("ATOM", "HETATM")):

            try:

                atom_name = line[12:16].strip()

                element = line[77:79].strip().upper()

                if not element:

                    element = re.sub(r"[^A-Za-z]", "", atom_name[:2]).upper()[:1]



                x = float(line[30:38])

                y = float(line[38:46])

                z = float(line[46:54])

                atoms.append(

                    {

                        "name": atom_name,

                        "element": element,

                        "coord": np.array([x, y, z], dtype=float),

                    }

                )

            except (ValueError, IndexError):

                continue

    return atoms





def parse_pdb_ligand_atoms(pdb_file, ligand_name, chain_id=None):

    parser = PDBParser(QUIET=True)

    structure = parser.get_structure("ligand", str(pdb_file))

    atoms = []



    for model in structure:

        for chain in model:

            if chain_id and chain.id != chain_id:

                continue

            for residue in chain:

                if residue.resname.strip() != ligand_name.strip():

                    continue

                for atom in residue:

                    element = (atom.element or "").strip().upper()

                    if not element:

                        element = re.sub(r"[^A-Za-z]", "", atom.name).upper()[:1]

                    atoms.append(

                        {

                            "name": atom.name.strip(),

                            "element": element,

                            "coord": np.array(atom.coord, dtype=float),

                        }

                    )

    return atoms





def kabsch_rmsd(reference_coords, predicted_coords):

    reference = np.asarray(reference_coords, dtype=float)

    predicted = np.asarray(predicted_coords, dtype=float)



    if reference.shape != predicted.shape or len(reference) < 3:

        raise ValueError("At least 3 matching atoms are required for RMSD.")



    ref_center = reference.mean(axis=0)

    pred_center = predicted.mean(axis=0)



    ref = reference - ref_center

    pred = predicted - pred_center



    covariance = pred.T @ ref

    v, s, wt = np.linalg.svd(covariance)

    d = np.sign(np.linalg.det(v @ wt))

    rotation = v @ np.diag([1.0, 1.0, d]) @ wt



    aligned = pred @ rotation

    diff = ref - aligned

    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))





def calculate_redocking_rmsd(

    reference_ligand_pdb,

    docked_pdbqt,

):

    """

    Calculate RMSD between the experimental reference ligand and

    the best redocked ligand pose.



    The experimental ligand is converted to an RDKit molecule and

    atom correspondence is determined using molecular structure

    rather than relying on PDB atom names.

    """



    print("\nCalculating redocking RMSD...")



    # ---------------------------------------------------------

    # 1. Load experimental ligand

    # ---------------------------------------------------------

    reference_mol = Chem.MolFromPDBFile(

        str(reference_ligand_pdb),

        removeHs=False,

        sanitize=False

    )



    if reference_mol is None:

        raise ValueError(

            "Could not load the experimental reference ligand with RDKit."

        )



    # Try sanitization after loading.

    try:

        Chem.SanitizeMol(reference_mol)

    except Exception:

        # Some PDB ligands do not contain enough connectivity

        # information for complete sanitization.

        print(

            "Warning: experimental ligand could not be fully sanitized."

        )



    # ---------------------------------------------------------

    # 2. Read docked PDBQT coordinates

    # ---------------------------------------------------------

    docked_atoms = parse_pdbqt_atoms(docked_pdbqt)



    if not docked_atoms:

        raise ValueError(

            "No atoms were found in the docked PDBQT file."

        )



    # Use MODEL 1 only.

    print(

        f"Experimental ligand atoms: {reference_mol.GetNumAtoms()}"

    )



    print(

        f"Docked ligand atoms: {len(docked_atoms)}"

    )



    # ---------------------------------------------------------

    # 3. Extract experimental coordinates

    # ---------------------------------------------------------

    conformer = reference_mol.GetConformer()



    reference_coordinates = []



    for atom_index in range(reference_mol.GetNumAtoms()):

        position = conformer.GetAtomPosition(atom_index)



        reference_coordinates.append(

            np.array(

                [position.x, position.y, position.z],

                dtype=float

            )

        )



    # ---------------------------------------------------------

    # 4. Match atoms using element + spatial correspondence

    #

    # First attempt:

    # match atoms by element and nearest coordinates.

    #

    # This does NOT require PDB atom names to match.

    # ---------------------------------------------------------

    from scipy.optimize import linear_sum_assignment



    reference_elements = [

        atom.GetSymbol().upper()

        for atom in reference_mol.GetAtoms()

    ]



    docked_elements = [

        atom["element"].upper()

        for atom in docked_atoms

    ]



    # Group atoms by element.

    reference_by_element = {}



    for i, element in enumerate(reference_elements):

        reference_by_element.setdefault(element, []).append(i)



    docked_by_element = {}



    for i, element in enumerate(docked_elements):

        docked_by_element.setdefault(element, []).append(i)



    matched_reference = []

    matched_docked = []



    # ---------------------------------------------------------

    # Match each element type independently.

    # ---------------------------------------------------------

    for element in reference_by_element:



        if element not in docked_by_element:

            continue



        ref_indices = reference_by_element[element]

        dock_indices = docked_by_element[element]



        ref_coords = np.array(

            [reference_coordinates[i] for i in ref_indices],

            dtype=float

        )



        dock_coords = np.array(

            [

                docked_atoms[i]["coord"]

                for i in dock_indices

            ],

            dtype=float

        )



        if len(ref_indices) == 0 or len(dock_indices) == 0:

            continue



        # Distance matrix.

        distances = np.linalg.norm(

            ref_coords[:, None, :] -

            dock_coords[None, :, :],

            axis=2

        )



        rows, cols = linear_sum_assignment(distances)



        for row, col in zip(rows, cols):

            matched_reference.append(

                ref_indices[row]

            )



            matched_docked.append(

                dock_indices[col]

            )



    # ---------------------------------------------------------

    # 5. Verify enough atoms were matched.

    # ---------------------------------------------------------

    matched_count = len(matched_reference)



    print(

        f"Matched ligand atoms: {matched_count}"

    )



    if matched_count < 3:



        # Provide useful diagnostics instead of a vague error.

        print("\nReference ligand elements:")

        print(reference_elements)



        print("\nDocked ligand elements:")

        print(docked_elements)



        raise ValueError(

            "Could not obtain at least 3 matching ligand atoms "

            "for RMSD. The experimental and docked ligand "

            "atom representations are incompatible."

        )



    # ---------------------------------------------------------

    # 6. Build coordinate arrays

    # ---------------------------------------------------------

    reference_coords = np.array(

        [

            reference_coordinates[i]

            for i in matched_reference

        ],

        dtype=float

    )



    docked_coords = np.array(

        [

            docked_atoms[i]["coord"]

            for i in matched_docked

        ],

        dtype=float

    )



    # ---------------------------------------------------------

    # 7. Calculate Kabsch RMSD

    # ---------------------------------------------------------

    rmsd = kabsch_rmsd(

        reference_coords,

        docked_coords

    )



    print(

        f"Redocking RMSD: {rmsd:.3f} Å"

    )



    return (

        float(rmsd),

        matched_count,

        "element-based correspondence + Kabsch"

    )





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

    run_dir = RESULTS_DIR / safe_name(protein_name) / safe_name(compound_name)

    run_dir.mkdir(parents=True, exist_ok=True)



    clean_pdb = run_dir / "protein_clean.pdb"

    sdf_file = run_dir / "ligand_3d.sdf"

    ligand_pdbqt = run_dir / "ligand.pdbqt"

    receptor_prefix = run_dir / "receptor"

    output_pdbqt = run_dir / "docking_output.pdbqt"



    selected_chain, _ = extract_and_clean_protein(

        pdb_file, clean_pdb, chain_id

    )



    smiles_to_sdf(smiles, sdf_file)

    prepare_ligand(sdf_file, ligand_pdbqt)

    receptor = prepare_receptor(clean_pdb, receptor_prefix)



    if site_mode == "reference":

        box = calculate_reference_ligand_box(

            pdb_file, reference_ligand, padding=8.0, chain_id=selected_chain

        )

    else:

        box = calculate_blind_box(clean_pdb)



    vina_output = run_vina(

        receptor,

        ligand_pdbqt,

        box,

        output_pdbqt,

        exhaustiveness,

        num_modes,

    )



    score = parse_vina_score(output_pdbqt)



    return {

        "protein": protein_name,

        "compound": compound_name,

        "binding_affinity": score,

        "selected_chain": selected_chain,

        "output": str(output_pdbqt),

        "protein_pdb": str(clean_pdb),

        "ligand_pdbqt": str(ligand_pdbqt),

        "box": box,

        "vina_output": vina_output,

    }





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

        / safe_name(protein_name)

        / f"redocking_{safe_name(reference_ligand)}"

    )

    run_dir.mkdir(parents=True, exist_ok=True)



    clean_pdb = run_dir / "protein_clean.pdb"

    reference_ligand_pdb = run_dir / "reference_ligand_experimental.pdb"

    reference_ligand_pdbqt = run_dir / "reference_ligand.pdbqt"

    receptor_prefix = run_dir / "receptor"

    output_pdbqt = run_dir / "redocking_output.pdbqt"



    selected_chain, _ = extract_and_clean_protein(

        pdb_file, clean_pdb, chain_id

    )



    extract_reference_ligand(

        pdb_file,

        reference_ligand_pdb,

        reference_ligand,

        selected_chain,

    )



    reference_mol = load_reference_ligand_as_rdkit(reference_ligand_pdb)

    prepare_rdkit_ligand(reference_mol, reference_ligand_pdbqt)



    receptor = prepare_receptor(clean_pdb, receptor_prefix)



    box = calculate_reference_ligand_box(

        pdb_file,

        reference_ligand,

        padding=8.0,

        chain_id=selected_chain,

    )



    vina_output = run_vina(

        receptor,

        reference_ligand_pdbqt,

        box,

        output_pdbqt,

        exhaustiveness,

        num_modes,

    )



    score = parse_vina_score(output_pdbqt)



    rmsd, matched_atoms, rmsd_method = calculate_redocking_rmsd(
        reference_ligand_pdb,
        output_pdbqt,
    )



    return {

        "protein": protein_name,

        "reference_ligand": reference_ligand,

        "binding_affinity": score,

        "rmsd": rmsd,

        "matched_atoms": matched_atoms,

        "rmsd_method": rmsd_method,

        "selected_chain": selected_chain,

        "experimental_ligand": str(reference_ligand_pdb),

        "predicted_output": str(output_pdbqt),

        "protein_pdb": str(clean_pdb),

        "receptor_pdbqt": str(receptor),

        "reference_ligand_pdbqt": str(reference_ligand_pdbqt),

        "box": box,

        "vina_output": vina_output,

    }





st.set_page_config(

    page_title="DistDockNet",

    page_icon="🧬",

    layout="wide",

)



st.title("🧬 DistDockNet")

st.caption("Protein–compound molecular docking and redocking validation")
st.caption("Protein structures are permanently stored in SQLite (distdocknet.db). Working PDB files are created only when the docking tools need them.")



if not VINA_EXE.exists():

    st.error(

        f"vina.exe was not found at:\n{VINA_EXE}\n\n"

        "Place vina.exe in the same folder as app.py."

    )



st.sidebar.header("Docking settings")



mode = st.sidebar.selectbox(

    "Workflow",

    ["Standard docking", "Redocking validation"],

)



site_mode = "Blind docking"

if mode == "Standard docking":

    site_mode = st.sidebar.selectbox(

        "Docking mode",

        ["Blind docking", "Reference ligand"],

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



if mode == "Redocking validation":

    st.info(

        "Redocking validation uses the ligand already present in the experimental "

        "PDB structure. The ligand is extracted, the protein is cleaned, the same "

        "ligand is redocked, and the predicted pose is compared with the experimental pose."

    )



    st.subheader("Redocking validation input")



    c1, c2 = st.columns(2)

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

        st.write("Expected example: MMP-13 / 2OW9 / SP6")



    st.divider()

    if st.button(

        "🧪 Run redocking validation",

        type="primary",

        use_container_width=True,

    ):

        if not VINA_EXE.exists():

            st.error("vina.exe is missing.")

            st.stop()



        if not protein_name.strip() or not pdb_id.strip():

            st.error("Enter a protein name and PDB ID.")

            st.stop()



        if not reference_ligand.strip():

            st.error("Enter the reference ligand CCD code.")

            st.stop()



        try:

            redocking_run_dir = (
                RESULTS_DIR
                / safe_name(protein_name)
                / f"redocking_{safe_name(reference_ligand)}"
            )

            with st.spinner(f"Getting {pdb_id.upper()} from database/RCSB..."):
                pdb_file = fetch_pdb(
                    pdb_id,
                    redocking_run_dir,
                )



            with st.spinner(

                f"Redocking {reference_ligand.upper()} into {protein_name}..."

            ):

                result = run_redocking_validation(

                    protein_name=protein_name,

                    pdb_file=pdb_file,

                    reference_ligand=reference_ligand,

                    chain_id=chain_id.strip() or "A",

                    exhaustiveness=exhaustiveness,

                    num_modes=num_modes,

                )



            st.success("Redocking completed successfully.")



            c1, c2, c3 = st.columns(3)

            with c1:

                st.metric(

                    "Vina score",

                    f"{result['binding_affinity']:.3f} kcal/mol",

                )

            with c2:

                st.metric(

                    "Pose RMSD",

                    f"{result['rmsd']:.3f} Å",

                )

            with c3:

                st.metric(

                    "Matched atoms",

                    str(result["matched_atoms"]),

                )



            if result["rmsd"] <= 2.0:

                st.success(

                    "RMSD ≤ 2.0 Å: the redocked pose reproduces the experimental "

                    "ligand position reasonably well."

                )

            else:

                st.warning(

                    "RMSD > 2.0 Å: the redocked pose differs substantially from "

                    "the experimental pose. Check the docking box, ligand preparation, "

                    "protonation, and receptor preparation."

                )



            st.write("**RMSD calculation:**", result["rmsd_method"])

            st.write("**Experimental ligand:**", result["experimental_ligand"])

            st.write("**Predicted docking output:**", result["predicted_output"])

            st.write("**Docking box:**", result["box"])



            st.session_state["redocking_result"] = result



        except Exception as exc:

            st.error("Redocking validation failed.")

            st.exception(exc)



    if "redocking_result" in st.session_state:

        result = st.session_state["redocking_result"]



        st.divider()

        st.subheader("Redocking files")



        file_rows = pd.DataFrame(

            [

                {

                    "File": "Experimental ligand",

                    "Path": result["experimental_ligand"],

                },

                {

                    "File": "Prepared reference ligand",

                    "Path": result["reference_ligand_pdbqt"],

                },

                {

                    "File": "Clean protein",

                    "Path": result["protein_pdb"],

                },

                {

                    "File": "Receptor PDBQT",

                    "Path": result["receptor_pdbqt"],

                },

                {

                    "File": "Predicted pose",

                    "Path": result["predicted_output"],

                },

            ]

        )

        st.dataframe(file_rows, use_container_width=True, hide_index=True)



else:

    st.subheader("1. Proteins")
    st.caption("Enter a PDB ID. DistDockNet checks SQLite first and fetches from RCSB PDB only when the structure is not already stored.")



    protein_count = st.number_input(

        "Number of proteins",

        min_value=1,

        max_value=20,

        value=1,

        step=1,

    )



    proteins = []



    for i in range(int(protein_count)):

        st.markdown(f"**Protein {i + 1}**")

        c1, c2 = st.columns(2)



        with c1:
            name = st.text_input(
                "Protein name",
                value=f"Protein_{i + 1}",
                key=f"protein_name_{i}",
            )

        with c2:
            pdb_id = st.text_input(
                "PDB ID",
                placeholder="Example: 2OW9",
                key=f"pdb_id_{i}",
            )

        proteins.append(
            {
                "name": name,
                "pdb_id": pdb_id,
            }
        )

    st.subheader("2. Compounds")



    compound_count = st.number_input(

        "Number of compounds",

        min_value=1,

        max_value=20,

        value=1,

        step=1,

    )



    compounds = []



    for i in range(int(compound_count)):

        st.markdown(f"**Compound {i + 1}**")

        name = st.text_input(

            "Compound name",

            value=f"Compound_{i + 1}",

            key=f"compound_name_{i}",

        )

        smiles = st.text_area(

            "SMILES",

            placeholder="Paste the compound SMILES here",

            key=f"compound_smiles_{i}",

            height=70,

        )



        compounds.append(

            {

                "name": name,

                "smiles": smiles,

            }

        )



    reference_ligand = ""

    if site_mode == "Reference ligand":

        reference_ligand = st.text_input(

            "Reference ligand CCD code",

            placeholder="Example: SP6",

        )



    st.divider()

    st.subheader("3. Run docking")



    if st.button(

        "🚀 Start docking",

        type="primary",

        use_container_width=True,

    ):

        if not VINA_EXE.exists():

            st.error("vina.exe is missing.")

            st.stop()



        valid_compounds = [

            c for c in compounds

            if c["name"].strip() and c["smiles"].strip()

        ]



        if not valid_compounds:

            st.error("Enter at least one valid compound.")

            st.stop()



        resolved_proteins = []



        try:

            for protein in proteins:

                if not protein["name"].strip():

                    raise ValueError("Every protein needs a name.")



                if protein["pdb_id"].strip():
                    working_directory = (
                        RESULTS_DIR
                        / safe_name(protein["name"])
                        / "_protein_input"
                    )

                    with st.spinner(
                        f"Getting {protein['pdb_id'].upper()} from database/RCSB..."
                    ):
                        pdb_file = fetch_pdb(
                            protein["pdb_id"],
                            working_directory,
                        )
                else:
                    raise ValueError(
                        f"Protein '{protein['name']}' needs a PDB ID."
                    )

                resolved_proteins.append(

                    {

                        "name": protein["name"],

                        "pdb_file": pdb_file,

                    }

                )



        except Exception as exc:

            st.error(str(exc))

            st.stop()



        jobs = list(product(resolved_proteins, valid_compounds))



        st.info(

            f"Running {len(jobs)} docking combination(s): "

            f"{len(resolved_proteins)} protein(s) × "

            f"{len(valid_compounds)} compound(s)."

        )



        progress = st.progress(0)

        status_box = st.empty()

        results = []



        for index, (protein, compound) in enumerate(jobs, start=1):

            status_box.write(

                f"Docking {index}/{len(jobs)}: "

                f"{protein['name']} + {compound['name']}"

            )



            try:

                result = run_standard_docking(

                    protein_name=protein["name"],

                    pdb_file=protein["pdb_file"],

                    compound_name=compound["name"],

                    smiles=compound["smiles"],

                    chain_id=chain_id.strip() or "A",

                    site_mode=(

                        "reference"

                        if site_mode == "Reference ligand"

                        else "blind"

                    ),

                    reference_ligand=reference_ligand,

                    exhaustiveness=exhaustiveness,

                    num_modes=num_modes,

                )



                results.append(

                    {

                        "Protein": result["protein"],

                        "Compound": result["compound"],

                        "Binding affinity (kcal/mol)": result[

                            "binding_affinity"

                        ],

                        "Chain": result["selected_chain"],

                        "Docking output": result["output"],

                    }

                )



            except Exception as exc:

                results.append(

                    {

                        "Protein": protein["name"],

                        "Compound": compound["name"],

                        "Binding affinity (kcal/mol)": None,

                        "Chain": None,

                        "Docking output": "",

                        "Error": str(exc),

                    }

                )



            progress.progress(index / len(jobs))



        status_box.success("Docking batch completed.")



        results_df = pd.DataFrame(results)

        output_csv = RESULTS_DIR / "ui_results.csv"

        results_df.to_csv(output_csv, index=False)

        st.session_state["results_df"] = results_df



        st.success(

            f"Completed {len(jobs)} docking combination(s)."

        )



    if "results_df" in st.session_state:

        st.divider()

        st.subheader("4. Docking results")



        results_df = st.session_state["results_df"]



        successful = results_df[

            results_df["Binding affinity (kcal/mol)"].notna()

        ].copy()



        if not successful.empty:

            successful = successful.sort_values(

                "Binding affinity (kcal/mol)"

            )



            st.dataframe(

                successful,

                use_container_width=True,

                hide_index=True,

            )



            st.subheader("Best docking score")

            best = successful.iloc[0]



            st.metric(

                label=(

                    f"{best['Protein']} + "

                    f"{best['Compound']}"

                ),

                value=(

                    f"{best['Binding affinity (kcal/mol)']:.3f} "

                    "kcal/mol"

                ),

            )



        failed = results_df[

            results_df["Binding affinity (kcal/mol)"].isna()

        ]



        if not failed.empty:

            st.warning(

                f"{len(failed)} docking job(s) failed."

            )

            st.dataframe(

                failed,

                use_container_width=True,

                hide_index=True,

            )



        csv_data = results_df.to_csv(index=False).encode("utf-8")



        st.download_button(

            "Download results CSV",

            data=csv_data,

            file_name="distdocknet_ui_results.csv",

            mime="text/csv",

        )



st.divider()

st.caption(

    "Scores are Vina docking predictions and should not be interpreted "

    "as experimental binding affinities or pharmacological potency."

)
