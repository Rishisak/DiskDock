"""
docking_bridge.py - Integrates DistDockNet's Molecular Docking Engine
into the Disease Compound Analysis Website and API.
"""

import json
import logging
import os
import re
import sys
import threading
import time
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

# Add project roots to path so DistDockNet modules can be imported
CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent if CURRENT_DIR.name == "Disease-Compound-Analysis" else CURRENT_DIR
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

logger = logging.getLogger(__name__)

# Paths
VINA_EXE = ROOT_DIR / "vina.exe"
PROTEIN_CACHE = ROOT_DIR / "protein" / "cache"
RESULTS_DIR = ROOT_DIR / "results"
UPLOAD_DIR = ROOT_DIR / "results" / "uploads"

PROTEIN_CACHE.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# Curated High-Affinity / Co-Crystallized Structure PDB Mapping
CURATED_TARGET_PDBS = {
    "PPARG": "2P4Y",
    "TNF": "2AZ5",
    "EGFR": "1M17",
    "BRAF": "4MNE",
    "CDK1": "6GU2",
    "CDK2": "1HCK",
    "CDK4": "2W96",
    "CDK6": "1BLX",
    "MMP13": "2OW9",
    "MMP9": "1L6J",
    "MMP2": "1CK7",
    "AKT1": "3O96",
    "AKT2": "1MRV",
    "MAPK1": "6QGN",
    "MAPK3": "4QTB",
    "MAPK14": "1A9U",
    "PIK3CA": "4JPS",
    "MTOR": "4JSP",
    "INSR": "1IR3",
    "INS": "4INS",
    "VEGFA": "1VPF",
    "KDR": "2XIR",
    "ESR1": "1ERE",
    "AR": "1E3G",
    "ACE2": "1R42",
    "PTGS2": "5KIR",
    "COX2": "5KIR",
    "PTGS1": "3N8Z",
    "COX1": "3N8Z",
    "HMGCR": "1HW8",
    "DPP4": "4A5S",
    "GLP1R": "5VEW",
    "STAT3": "6TLC",
    "RELA": "1NFI",
    "NFKB1": "1NFI",
    "TP53": "1TUP",
    "GSK3B": "1H8F",
    "JAK2": "3HCT",
    "PARP1": "4UND",
    "BCL2": "4LVT",
    "HSP90AA1": "2VCI",
    "IL6": "1ALU",
    "IL1B": "9ILB",
    "NOS2": "1NOD",
    "NOS3": "1M9Q",
    "SIRT1": "4ZZJ",
    "PRKAA1": "4CFE",
    "KCNJ11": "7WIT",
    "ABCC8": "6PZB",
}

# In-memory docking job status tracking
_jobs_lock = threading.Lock()
_docking_jobs: Dict[str, Dict[str, Any]] = {}


def safe_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", (value or "").strip())
    return cleaned.strip("._") or "unnamed"


def resolve_target_pdb(target_symbol: str, target_name: str = "") -> Dict[str, str]:
    """Resolve a gene symbol or protein target to a 4-letter PDB ID."""
    sym = (target_symbol or "").strip().upper()
    if sym in CURATED_TARGET_PDBS:
        return {
            "target": sym,
            "pdb_id": CURATED_TARGET_PDBS[sym],
            "source": "Curated High-Res Crystal Structure",
        }

    try:
        from protein.lookup import search_pdb
        pdb_id = search_pdb(sym)
        if pdb_id and len(pdb_id) == 4:
            return {
                "target": sym,
                "pdb_id": pdb_id.upper(),
                "source": "RCSB PDB Dynamic Query",
            }
    except Exception as exc:
        logger.warning(f"RCSB search for {sym} failed: {exc}")

    if target_name and target_name.strip():
        try:
            from protein.lookup import search_pdb
            pdb_id = search_pdb(target_name.strip())
            if pdb_id and len(pdb_id) == 4:
                return {
                    "target": sym,
                    "pdb_id": pdb_id.upper(),
                    "source": "RCSB PDB Dynamic Query (Name)",
                }
        except Exception:
            pass

    return {
        "target": sym,
        "pdb_id": "",
        "source": "Unresolved",
    }


def fetch_pdb_file(pdb_id: str) -> Path:
    """Download or return cached PDB file from RCSB PDB."""
    clean_id = pdb_id.strip().upper()
    if not re.fullmatch(r"[0-9A-Za-z]{4}", clean_id):
        raise ValueError(f"Invalid PDB ID format: '{pdb_id}'. Must be exactly 4 alphanumeric characters.")

    target_path = PROTEIN_CACHE / f"{clean_id}.pdb"
    if target_path.exists() and target_path.stat().st_size > 500:
        return target_path

    url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    try:
        urllib.request.urlretrieve(url, target_path)
    except Exception as exc:
        if target_path.exists():
            target_path.unlink()
        raise RuntimeError(f"Could not retrieve PDB structure {clean_id} from RCSB: {exc}")

    return target_path


def pdbqt_pose_to_pdb(pdbqt_text: str) -> str:
    """
    Convert Vina PDBQT atom coordinates into clean standard PDB HETATM lines.
    Allows 3Dmol.js to render the docked ligand without PDBQT parser issues.
    """
    lines = []
    serial = 1

    for line in pdbqt_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            atom_name = line[12:16].strip() or "X"
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            # Parse element
            raw_element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not raw_element:
                alpha_chars = "".join(c for c in atom_name if c.isalpha())
                raw_element = alpha_chars[:2].upper() if len(alpha_chars) >= 2 else (alpha_chars[:1].upper() or "C")

            if raw_element == "CL":
                elem = "Cl"
            elif raw_element == "BR":
                elem = "Br"
            else:
                elem = raw_element[:1].upper()

            pdb_line = (
                f"HETATM"
                f"{serial:5d} "
                f"{atom_name:<4s}"
                f" LIG A   1    "
                f"{x:8.3f}"
                f"{y:8.3f}"
                f"{z:8.3f}"
                f"  1.00  0.00          "
                f"{elem:>2s}"
            )
            lines.append(pdb_line)
            serial += 1
        except Exception:
            continue

    lines.append("END")
    return "\n".join(lines)


def extract_vina_poses(pdbqt_text: str) -> List[Dict[str, Any]]:
    """Parse AutoDock Vina multi-pose output and extract scores and PDB structures."""
    poses = []
    models = re.findall(r"MODEL\s+(\d+)(.*?)(?:ENDMDL|$)", pdbqt_text, re.DOTALL)

    for model_num, model_body in models:
        score_match = re.search(
            r"REMARK VINA RESULT:\s*([-+]?\d+(?:\.\d+)?)\s*([-+]?\d+(?:\.\d+)?)\s*([-+]?\d+(?:\.\d+)?)",
            model_body,
        )
        score = float(score_match.group(1)) if score_match else None
        rmsd_lb = float(score_match.group(2)) if score_match else 0.0
        rmsd_ub = float(score_match.group(3)) if score_match else 0.0

        model_pdbqt = f"MODEL {model_num}\n" + model_body.strip() + "\nENDMDL"
        model_pdb = pdbqt_pose_to_pdb(model_body)

        poses.append({
            "pose": int(model_num),
            "score": score,
            "rmsd_lb": rmsd_lb,
            "rmsd_ub": rmsd_ub,
            "pdbqt": model_pdbqt,
            "pdb": model_pdb,
        })

    return poses


def run_docking_pipeline(
    protein_name: str,
    pdb_id: str = "",
    uploaded_pdb_content: Optional[str] = None,
    compound_name: str = "Compound",
    smiles: str = "",
    chain_id: str = "A",
    exhaustiveness: int = 16,
    num_modes: int = 9,
    site_mode: str = "blind",
    reference_ligand: str = "",
    job_id: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Executes actual AutoDock Vina molecular docking using DistDockNet's modules.
    Returns complete results including binding affinities and 3D PDB coordinates.
    """
    if not VINA_EXE.exists():
        raise FileNotFoundError(f"AutoDock Vina executable not found at: {VINA_EXE}")

    if not smiles or not smiles.strip():
        raise ValueError("Valid chemical SMILES string is required for docking.")

    p_clean = safe_filename(protein_name)
    c_clean = safe_filename(compound_name)
    run_dir = RESULTS_DIR / p_clean / c_clean
    run_dir.mkdir(parents=True, exist_ok=True)

    # 1. Obtain receptor PDB
    if uploaded_pdb_content and len(uploaded_pdb_content.strip()) > 100:
        pdb_file = run_dir / "uploaded_input.pdb"
        pdb_file.write_text(uploaded_pdb_content, encoding="utf-8")
    elif pdb_id and pdb_id.strip():
        pdb_file = fetch_pdb_file(pdb_id.strip())
    else:
        # Try resolving by target name
        resolved = resolve_target_pdb(protein_name)
        if resolved.get("pdb_id"):
            pdb_id = resolved["pdb_id"]
            pdb_file = fetch_pdb_file(pdb_id)
        else:
            raise ValueError(f"No PDB ID provided or resolved for target '{protein_name}'.")

    # 2. Extract chain & clean protein
    from protein.structure import extract_chain
    from protein.preparation import clean_protein, prepare_receptor
    from compound.structure import smiles_to_3d
    from compound.preparation import prepare_ligand
    from docking.box import calculate_blind_docking_box, calculate_box_from_reference_ligand
    from docking.vina_runner import run_vina

    chain_pdb = run_dir / "protein_chain.pdb"
    clean_pdb = run_dir / "protein_clean.pdb"
    sdf_file = run_dir / "ligand_3d.sdf"
    ligand_pdbqt = run_dir / "ligand.pdbqt"
    receptor_prefix = run_dir / "receptor"
    output_pdbqt = run_dir / "docking_output.pdbqt"

    # Chain extraction
    try:
        chain_pdb_path, selected_chain_id = extract_chain(
            str(pdb_file),
            str(chain_pdb),
            chain_id.strip() or "A",
        )
    except Exception as exc:
        logger.warning(f"Chain extraction warning ({exc}), fallback to whole PDB")
        chain_pdb_path = pdb_file
        selected_chain_id = chain_id.strip() or "A"

    # Clean protein
    clean_protein(str(chain_pdb_path), str(clean_pdb))

    # Receptor PDBQT
    prepare_receptor(str(clean_pdb), str(receptor_prefix))
    receptor_pdbqt = receptor_prefix.with_suffix(".pdbqt")
    if not receptor_pdbqt.exists():
        # Fallback if Meeko created receptor.pdbqt without suffix
        potential = run_dir / "receptor.pdbqt"
        if potential.exists():
            receptor_pdbqt = potential

    # 3. Compound SMILES -> 3D SDF -> PDBQT
    smiles_to_3d(smiles.strip(), str(sdf_file))
    prepare_ligand(str(sdf_file), str(ligand_pdbqt))

    # 4. Docking Box
    if site_mode == "reference" and reference_ligand and reference_ligand.strip():
        box = calculate_box_from_reference_ligand(str(pdb_file), reference_ligand.strip())
    else:
        box = calculate_blind_docking_box(str(clean_pdb))

    # 5. Run AutoDock Vina
    vina_log = run_vina(
        vina_executable=str(VINA_EXE),
        receptor_pdbqt=str(receptor_pdbqt),
        ligand_pdbqt=str(ligand_pdbqt),
        output_pdbqt=str(output_pdbqt),
        box=box,
        exhaustiveness=int(exhaustiveness),
        num_modes=int(num_modes),
        seed=42,
    )

    # 6. Parse results & poses
    if not output_pdbqt.exists() or output_pdbqt.stat().st_size == 0:
        raise RuntimeError("AutoDock Vina finished without generating docking_output.pdbqt.")

    output_text = output_pdbqt.read_text(encoding="utf-8", errors="replace")
    poses = extract_vina_poses(output_text)
    best_affinity = poses[0]["score"] if (poses and poses[0]["score"] is not None) else None

    # Cleaned protein text for 3Dmol.js viewer
    protein_pdb_text = clean_pdb.read_text(encoding="utf-8", errors="replace") if clean_pdb.exists() else ""

    result_payload = {
        "status": "success",
        "job_id": job_id or f"{p_clean}_{c_clean}_{int(time.time())}",
        "protein_name": protein_name,
        "pdb_id": pdb_id.upper() if pdb_id else "",
        "compound_name": compound_name,
        "smiles": smiles,
        "chain_id": selected_chain_id,
        "binding_affinity": best_affinity,
        "best_affinity_kcal_mol": best_affinity,
        "num_modes": len(poses),
        "poses": poses,
        "box": box,
        "protein_pdb_text": protein_pdb_text,
        "output_pdbqt_path": str(output_pdbqt),
        "clean_pdb_path": str(clean_pdb),
        "ligand_pdbqt_path": str(ligand_pdbqt),
    }

    # 7. Persist to results CSV for Step 12 sync
    _save_to_csv_results(result_payload)

    return result_payload


def _save_to_csv_results(result: Dict[str, Any]):
    """Sync results into DistDockNet's batch_results.csv and ui_results.csv."""
    try:
        csv_row = {
            "protein": result.get("protein_name", ""),
            "Protein": result.get("protein_name", ""),
            "compound": result.get("compound_name", ""),
            "Compound": result.get("compound_name", ""),
            "affinity": result.get("binding_affinity"),
            "Binding affinity (kcal/mol)": result.get("binding_affinity"),
            "PDB ID": result.get("pdb_id", ""),
            "Chain": result.get("chain_id", "A"),
            "Docking output": result.get("output_pdbqt_path", ""),
            "Protein PDB": result.get("clean_pdb_path", ""),
            "Ligand PDBQT": result.get("ligand_pdbqt_path", ""),
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        }

        batch_csv = RESULTS_DIR / "batch_results.csv"
        ui_csv = RESULTS_DIR / "ui_results.csv"

        df_row = pd.DataFrame([csv_row])

        if batch_csv.exists() and batch_csv.stat().st_size > 0:
            df_row.to_csv(batch_csv, mode="a", header=False, index=False)
        else:
            df_row.to_csv(batch_csv, index=False)

        if ui_csv.exists() and ui_csv.stat().st_size > 0:
            df_row.to_csv(ui_csv, mode="a", header=False, index=False)
        else:
            df_row.to_csv(ui_csv, index=False)
    except Exception as exc:
        logger.warning(f"Could not persist result to CSV: {exc}")


def get_docking_history(limit: int = 50) -> List[Dict[str, Any]]:
    """Retrieve saved docking runs from batch_results.csv."""
    batch_csv = RESULTS_DIR / "batch_results.csv"
    if not batch_csv.exists():
        return []
    try:
        df = pd.read_csv(batch_csv)
        if df.empty:
            return []
        # Return most recent runs
        records = df.tail(limit).to_dict(orient="records")
        records.reverse()
        return records
    except Exception as exc:
        logger.warning(f"Error reading docking history: {exc}")
        return []

