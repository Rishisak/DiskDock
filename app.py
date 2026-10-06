
import io
import re
import subprocess
from pathlib import Path
from itertools import product

import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser, PDBIO, Select


# ============================================================
# DistDockNet - Basic Docking UI
# ============================================================

BASE_DIR = Path(__file__).resolve().parent
VINA_EXE = BASE_DIR / "vina.exe"
RESULTS_DIR = BASE_DIR / "results" / "ui_runs"
PROTEIN_CACHE = BASE_DIR / "protein" / "cache"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PROTEIN_CACHE.mkdir(parents=True, exist_ok=True)


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


def fetch_pdb(pdb_id):
    import urllib.request

    pdb_id = pdb_id.strip().upper()
    if not re.fullmatch(r"[0-9A-Za-z]{4}", pdb_id):
        raise ValueError("PDB ID must contain exactly 4 letters/numbers.")

    output = PROTEIN_CACHE / f"{pdb_id}.pdb"

    if output.exists() and output.stat().st_size > 0:
        return output

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        urllib.request.urlretrieve(url, output)
    except Exception as exc:
        if output.exists():
            output.unlink()
        raise RuntimeError(f"Could not download PDB {pdb_id}: {exc}") from exc

    return output


def get_largest_protein_chain(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_file))

    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("No model found in the protein structure.")

    counts = {}

    for chain in model:
        count = 0
        for residue in chain:
            if residue.resname.strip() in STANDARD_AMINO_ACIDS:
                count += 1
        if count > 0:
            counts[chain.id] = count

    if not counts:
        raise ValueError("No standard amino-acid protein chain was found.")

    return max(counts, key=counts.get)


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

    if chain_id not in available:
        chain_id = max(available, key=available.get)

    io_obj = PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(str(output_pdb), ProteinSelect(chain_id))

    return chain_id, available


def smiles_to_sdf(smiles, output_sdf):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    mol = Chem.AddHs(mol)

    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status != 0:
        raise RuntimeError("RDKit could not generate a 3D conformer.")

    AllChem.UFFOptimizeMolecule(mol)

    writer = Chem.SDWriter(str(output_sdf))
    writer.write(mol)
    writer.close()

    return output_sdf


def prepare_ligand(sdf_file, output_pdbqt):
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
    mol = supplier[0]

    if mol is None:
        raise ValueError("Could not read generated ligand structure.")

    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)

    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setups[0])

    if not is_ok:
        raise RuntimeError(f"Meeko ligand preparation failed: {error_msg}")

    Path(output_pdbqt).write_text(pdbqt_string, encoding="utf-8")
    return output_pdbqt


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

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

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

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    padding = 5.0

    return {
        "center_x": (min_x + max_x) / 2,
        "center_y": (min_y + max_y) / 2,
        "center_z": (min_z + max_z) / 2,
        "size_x": (max_x - min_x) + 2 * padding,
        "size_y": (max_y - min_y) + 2 * padding,
        "size_z": (max_z - min_z) + 2 * padding,
    }


def calculate_reference_ligand_box(pdb_file, ligand_resname, padding=8.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_file))

    coordinates = []

    for model in structure:
        for chain in model:
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


def run_vina(
    receptor,
    ligand,
    box,
    output_file,
    exhaustiveness=16,
    num_modes=10,
):
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

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        details = result.stderr or result.stdout or "Unknown Vina error."
        raise RuntimeError("Vina docking failed.\n\n" + details)

    return result.stdout


def parse_vina_score(pdbqt_file):
    pattern = re.compile(
        r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)"
    )

    for line in Path(pdbqt_file).read_text(
        encoding="utf-8",
        errors="ignore"
    ).splitlines():

        match = pattern.search(line)

        if match:
            return float(match.group(1))

    raise ValueError("Could not find VINA RESULT in docking output.")


def dock_one(
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
        / safe_name(protein_name)
        / safe_name(compound_name)
    )

    run_dir.mkdir(parents=True, exist_ok=True)

    clean_pdb = run_dir / "protein_clean.pdb"
    sdf_file = run_dir / "ligand_3d.sdf"
    ligand_pdbqt = run_dir / "ligand.pdbqt"
    receptor_prefix = run_dir / "receptor"
    output_pdbqt = run_dir / "docking_output.pdbqt"

    selected_chain, available = extract_and_clean_protein(
        pdb_file,
        clean_pdb,
        chain_id,
    )

    smiles_to_sdf(smiles, sdf_file)
    prepare_ligand(sdf_file, ligand_pdbqt)
    receptor = prepare_receptor(clean_pdb, receptor_prefix)

    if site_mode == "reference":
        box = calculate_reference_ligand_box(
            pdb_file,
            reference_ligand,
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


# ============================================================
# UI
# ============================================================

st.set_page_config(
    page_title="DistDockNet",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 DistDockNet")
st.caption(
    "Simple protein–compound molecular docking interface"
)

if not VINA_EXE.exists():
    st.error(
        f"vina.exe was not found at:\n{VINA_EXE}\n\n"
        "Place vina.exe in the same folder as app.py."
    )

st.sidebar.header("Docking settings")

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

st.subheader("1. Proteins")

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

        pdb_id = st.text_input(
            "PDB ID (optional)",
            placeholder="Example: 2OW9",
            key=f"pdb_id_{i}",
        )

    with c2:
        uploaded = st.file_uploader(
            "Or upload a PDB file",
            type=["pdb", "ent"],
            key=f"protein_file_{i}",
        )

    proteins.append(
        {
            "name": name,
            "pdb_id": pdb_id,
            "uploaded": uploaded,
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
        "Reference ligand residue name",
        placeholder="Example: ATP or other PDB ligand code",
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

            if protein["uploaded"] is not None:
                destination = (
                    PROTEIN_CACHE
                    / safe_name(protein["name"])
                    / "input.pdb"
                )
                destination.parent.mkdir(
                    parents=True,
                    exist_ok=True,
                )
                destination.write_bytes(
                    protein["uploaded"].getvalue()
                )
                pdb_file = destination

            elif protein["pdb_id"].strip():
                with st.spinner(
                    f"Getting {protein['pdb_id'].upper()}..."
                ):
                    pdb_file = fetch_pdb(
                        protein["pdb_id"]
                    )

            else:
                raise ValueError(
                    f"Protein '{protein['name']}' needs "
                    "either a PDB ID or a PDB file."
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

    jobs = list(
        product(
            resolved_proteins,
            valid_compounds,
        )
    )

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

            result = dock_one(
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
