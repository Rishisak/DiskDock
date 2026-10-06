"""Run and recover DistDockNet docking jobs.

This version:

1. Reads all 80 jobs from the batch CSV.
2. Preserves successful previous results.
3. Recovers completed Vina dockings from docking_output.pdbqt.
4. Retries only jobs that are genuinely unfinished.
5. Saves results after every job.
"""

import csv
from pathlib import Path
import re

from pipeline.run_pipeline import run_pipeline
from protein.remote import fetch_structure


PROJECT_ROOT = Path(__file__).resolve().parent

CSV_PATH = (
    PROJECT_ROOT /
    "DiscDockNet_80_batch.csv"
)

RESULTS_DIR = (
    PROJECT_ROOT /
    "results"
)

BATCH_RESULTS_FILE = (
    RESULTS_DIR /
    "batch_results.csv"
)

RESULT_MATRIX_FILE = (
    RESULTS_DIR /
    "binding_affinity_matrix.csv"
)

REQUIRED_COLUMNS = {
    "protein_name",
    "protein_pdb",
    "compound_name",
    "smiles",
}


# ======================================================
# Structure reference
# ======================================================

def parse_structure_reference(protein_pdb):

    reference = (
        protein_pdb
        .strip()
        .replace("\\", "/")
    )

    rcsb_match = re.fullmatch(
        r"protein/([A-Za-z0-9]{4})\.pdb",
        reference
    )

    alphafold_match = re.fullmatch(
        r"alphafold_([A-Za-z0-9]+)",
        reference,
        re.IGNORECASE
    )

    if rcsb_match:

        return (
            "RCSB",
            rcsb_match.group(1).upper()
        )

    if alphafold_match:

        return (
            "AlphaFold",
            alphafold_match.group(1).upper()
        )

    raise ValueError(
        f"Unknown protein structure format: "
        f"{protein_pdb!r}"
    )


# ======================================================
# Load original 80 jobs
# ======================================================

def load_batch_rows():

    if not CSV_PATH.exists():

        raise FileNotFoundError(
            f"Batch CSV not found: {CSV_PATH}"
        )

    with CSV_PATH.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as file:

        rows = list(
            csv.DictReader(file)
        )

    if not rows:

        raise ValueError(
            "The batch CSV contains no docking jobs."
        )

    missing_columns = (
        REQUIRED_COLUMNS -
        set(rows[0])
    )

    if missing_columns:

        raise ValueError(
            "CSV is missing required columns: "
            + ", ".join(
                sorted(missing_columns)
            )
        )

    return rows


# ======================================================
# Load current CSV results
# ======================================================

def load_previous_results():

    if not BATCH_RESULTS_FILE.exists():

        return {}

    with BATCH_RESULTS_FILE.open(
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as file:

        rows = list(
            csv.DictReader(file)
        )

    previous = {}

    for row in rows:

        key = (
            row.get(
                "protein",
                ""
            ).strip(),

            row.get(
                "compound",
                ""
            ).strip()
        )

        if key[0] and key[1]:

            previous[key] = row

    return previous


# ======================================================
# Extract Vina score from PDBQT
# ======================================================

def extract_vina_score(docking_output):

    docking_output = Path(
        docking_output
    )

    if not docking_output.exists():

        return None

    try:

        text = docking_output.read_text(
            encoding="utf-8",
            errors="ignore"
        )

    except Exception:

        return None

    patterns = [
        r"REMARK VINA RESULT:\s*"
        r"(-?\d+(?:\.\d+)?)",

        r"REMARK VINA RESULT\s+"
        r"(-?\d+(?:\.\d+)?)",
    ]

    for pattern in patterns:

        match = re.search(
            pattern,
            text
        )

        if match:

            try:

                return float(
                    match.group(1)
                )

            except ValueError:

                return None

    return None


# ======================================================
# Recover completed docking from disk
# ======================================================

def recover_completed_job(
    protein_name,
    compound_name
):

    work_dir = (
        RESULTS_DIR /
        safe_name(protein_name) /
        safe_name(compound_name)
    )

    docking_output = (
        work_dir /
        "docking_output.pdbqt"
    )

    score = extract_vina_score(
        docking_output
    )

    if score is None:

        return None

    metadata_file = (
        work_dir /
        "run_metadata.json"
    )

    selected_chain_id = ""

    if metadata_file.exists():

        try:

            import json

            metadata = json.loads(
                metadata_file.read_text(
                    encoding="utf-8"
                )
            )

            selected_chain_id = str(
                metadata.get(
                    "selected_chain_id",
                    ""
                )
            )

        except Exception:

            selected_chain_id = ""

    return {

        "protein": protein_name,

        "compound": compound_name,

        "binding_affinity": float(
            score
        ),

        "status": "SUCCESS",

        "selected_chain_id": (
            selected_chain_id
        ),

        "work_directory": str(
            work_dir
        ),

        "message": str(
            docking_output
        ),
    }


# ======================================================
# Safe directory names
# ======================================================

def safe_name(value):

    cleaned = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        value.strip()
    )

    return (
        cleaned.strip("._")
        or "unnamed"
    )


# ======================================================
# Write result files
# ======================================================

def write_result_files(
    batch_results
):

    RESULTS_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    fieldnames = [
        "protein",
        "compound",
        "binding_affinity",
        "status",
        "selected_chain_id",
        "work_directory",
        "message",
    ]

    # --------------------------------------------------
    # Remove duplicate protein-compound rows
    # Keep the latest result.
    # --------------------------------------------------

    unique_results = {}

    for result in batch_results:

        key = (
            result["protein"],
            result["compound"]
        )

        unique_results[key] = result

    ordered_results = []

    for key in unique_results:

        ordered_results.append(
            unique_results[key]
        )

    # --------------------------------------------------
    # Write batch_results.csv
    # --------------------------------------------------

    with BATCH_RESULTS_FILE.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames
        )

        writer.writeheader()

        writer.writerows(
            ordered_results
        )

    # --------------------------------------------------
    # Create matrix
    # --------------------------------------------------

    proteins = list(
        dict.fromkeys(
            result["protein"]
            for result in ordered_results
        )
    )

    compounds = list(
        dict.fromkeys(
            result["compound"]
            for result in ordered_results
        )
    )

    scores = {}

    for result in ordered_results:

        if (
            result["status"]
            .upper()
            == "SUCCESS"
        ):

            scores[
                (
                    result["protein"],
                    result["compound"]
                )
            ] = result[
                "binding_affinity"
            ]

    with RESULT_MATRIX_FILE.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=[
                "protein",
                *compounds
            ]
        )

        writer.writeheader()

        for protein in proteins:

            row = {
                "protein": protein
            }

            for compound in compounds:

                row[compound] = scores.get(
                    (
                        protein,
                        compound
                    ),
                    ""
                )

            writer.writerow(row)


# ======================================================
# Run batch
# ======================================================

def run_batch():

    print("\n========================================")
    print("       DistDockNet Smart Batch")
    print("========================================")

    rows = load_batch_rows()

    print(
        f"\nTotal jobs in CSV: {len(rows)}"
    )

    previous_results = (
        load_previous_results()
    )

    print(
        f"Existing CSV results: "
        f"{len(previous_results)}"
    )

    # --------------------------------------------------
    # Build complete result list
    # --------------------------------------------------

    results_by_key = {}

    # --------------------------------------------------
    # First: use existing successful CSV results
    # --------------------------------------------------

    for key, result in previous_results.items():

        if (
            result.get(
                "status",
                ""
            ).upper()
            == "SUCCESS"
        ):

            results_by_key[key] = {
                "protein": result[
                    "protein"
                ],

                "compound": result[
                    "compound"
                ],

                "binding_affinity": (
                    result.get(
                        "binding_affinity",
                        ""
                    )
                ),

                "status": "SUCCESS",

                "selected_chain_id": (
                    result.get(
                        "selected_chain_id",
                        ""
                    )
                ),

                "work_directory": (
                    result.get(
                        "work_directory",
                        ""
                    )
                ),

                "message": (
                    result.get(
                        "message",
                        ""
                    )
                ),
            }

    # --------------------------------------------------
    # Second: recover successful Vina outputs
    # --------------------------------------------------

    print(
        "\nChecking existing docking outputs..."
    )

    recovered_count = 0

    for row in rows:

        protein_name = (
            row["protein_name"].strip()
        )

        compound_name = (
            row["compound_name"].strip()
        )

        key = (
            protein_name,
            compound_name
        )

        # Already known success
        if key in results_by_key:

            continue

        recovered = recover_completed_job(
            protein_name,
            compound_name
        )

        if recovered is not None:

            results_by_key[key] = recovered

            recovered_count += 1

            print(
                "RECOVERED:",
                protein_name,
                "+",
                compound_name,
                "->",
                recovered[
                    "binding_affinity"
                ]
            )

    # --------------------------------------------------
    # Count recovered/successful jobs
    # --------------------------------------------------

    successful_count = sum(
        result["status"].upper()
        == "SUCCESS"

        for result in
        results_by_key.values()
    )

    print(
        f"\nSuccessful jobs preserved/recovered: "
        f"{successful_count}"
    )

    print(
        f"Newly recovered from Vina outputs: "
        f"{recovered_count}"
    )

    # --------------------------------------------------
    # Determine jobs that really need to run
    # --------------------------------------------------

    jobs_to_run = []

    for row in rows:

        protein_name = (
            row["protein_name"].strip()
        )

        compound_name = (
            row["compound_name"].strip()
        )

        key = (
            protein_name,
            compound_name
        )

        existing = results_by_key.get(
            key
        )

        if existing is not None:

            if (
                existing["status"]
                .upper()
                == "SUCCESS"
            ):

                print(
                    "SKIP SUCCESS:",
                    protein_name,
                    "+",
                    compound_name
                )

                continue

        jobs_to_run.append(row)

    # --------------------------------------------------
    # Save recovered results before new jobs
    # --------------------------------------------------

    write_result_files(
        list(results_by_key.values())
    )

    print("\n========================================")
    print("             JOB SELECTION")
    print("========================================")

    print(
        f"Total jobs      : {len(rows)}"
    )

    print(
        f"Successful      : {successful_count}"
    )

    print(
        f"Jobs to run     : {len(jobs_to_run)}"
    )

    # --------------------------------------------------
    # Nothing left
    # --------------------------------------------------

    if not jobs_to_run:

        print(
            "\nAll docking jobs are complete."
        )

        return

    # --------------------------------------------------
    # Run unfinished jobs
    # --------------------------------------------------

    for index, row in enumerate(
        jobs_to_run,
        start=1
    ):

        protein_name = (
            row["protein_name"].strip()
        )

        compound_name = (
            row["compound_name"].strip()
        )

        key = (
            protein_name,
            compound_name
        )

        print("\n" + "=" * 60)

        print(
            f"RUNNING JOB "
            f"{index}/{len(jobs_to_run)}"
        )

        print(
            f"{protein_name} + "
            f"{compound_name}"
        )

        print("=" * 60)

        try:

            source, structure_id = (
                parse_structure_reference(
                    row["protein_pdb"]
                )
            )

            print(
                f"\nStructure source: "
                f"{source}"
            )

            print(
                f"Structure ID: "
                f"{structure_id}"
            )

            # ------------------------------------------
            # Fetch protein
            # ------------------------------------------

            protein_file = fetch_structure(
                source,
                structure_id
            )

            # ------------------------------------------
            # Run complete pipeline
            # ------------------------------------------

            result = run_pipeline(

                protein_pdb=protein_file,

                compound_smiles=(
                    row["smiles"].strip()
                ),

                protein_name=protein_name,

                compound_name=compound_name,

                chain_id="A",

                site_mode="blind",

                results_dir=RESULTS_DIR,
            )

            # ------------------------------------------
            # Mark SUCCESS
            # ------------------------------------------

            results_by_key[key] = {

                "protein": protein_name,

                "compound": compound_name,

                "binding_affinity": float(
                    result[
                        "binding_affinity"
                    ]
                ),

                "status": "SUCCESS",

                "selected_chain_id": (
                    result[
                        "selected_chain_id"
                    ]
                ),

                "work_directory": (
                    result[
                        "work_directory"
                    ]
                ),

                "message": (
                    result[
                        "docking_output"
                    ]
                ),
            }

            print(
                "\nJOB COMPLETED SUCCESSFULLY"
            )

            print(
                "Binding affinity:",
                result[
                    "binding_affinity"
                ]
            )

        except Exception as error:

            # ------------------------------------------
            # IMPORTANT:
            # Check whether Vina actually produced
            # an output despite a later error.
            # ------------------------------------------

            recovered = recover_completed_job(
                protein_name,
                compound_name
            )

            if recovered is not None:

                results_by_key[key] = recovered

                print(
                    "\nDOCKING OUTPUT FOUND."
                )

                print(
                    "The job completed docking "
                    "even though a later pipeline "
                    "step raised an error."
                )

                print(
                    "Recovered affinity:",
                    recovered[
                        "binding_affinity"
                    ]
                )

            else:

                message = (
                    f"{type(error).__name__}: "
                    f"{error}"
                )

                results_by_key[key] = {

                    "protein": protein_name,

                    "compound": compound_name,

                    "binding_affinity": "",

                    "status": "FAILED",

                    "selected_chain_id": "",

                    "work_directory": "",

                    "message": message,
                }

                print(
                    "\nJOB FAILED:"
                )

                print(
                    message
                )

        # ----------------------------------------------
        # Save immediately after every job
        # ----------------------------------------------

        write_result_files(
            list(
                results_by_key.values()
            )
        )

    # --------------------------------------------------
    # Final summary
    # --------------------------------------------------

    final_results = list(
        results_by_key.values()
    )

    successful = sum(
        result["status"].upper()
        == "SUCCESS"

        for result in final_results
    )

    failed = sum(
        result["status"].upper()
        == "FAILED"

        for result in final_results
    )

    print("\n" + "=" * 60)

    print(
        "BATCH PROCESSING COMPLETED"
    )

    print("=" * 60)

    print(
        f"Total jobs : {len(rows)}"
    )

    print(
        f"Successful : {successful}"
    )

    print(
        f"Failed     : {failed}"
    )

    print(
        f"Results    : {BATCH_RESULTS_FILE}"
    )

    print(
        f"Matrix     : {RESULT_MATRIX_FILE}"
    )


if __name__ == "__main__":
    run_batch()
