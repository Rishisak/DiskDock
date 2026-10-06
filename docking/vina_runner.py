from pathlib import Path
import re
import subprocess


VINA_RESULT_PATTERN = re.compile(
    r"^REMARK VINA RESULT:\s*([-+]?\d+(?:\.\d+)?)", re.MULTILINE
)


def extract_best_affinity(output_pdbqt, vina_stdout=""):
    """Read Vina's best affinity from its PDBQT output (or stdout fallback)."""
    output_path = Path(output_pdbqt)
    text = output_path.read_text(encoding="utf-8", errors="replace")
    match = VINA_RESULT_PATTERN.search(text)

    if match is None:
        # Vina's console table is retained as a fallback for version differences.
        match = re.search(
            r"^\s*1\s+([-+]?\d+(?:\.\d+)?)\s+", vina_stdout, re.MULTILINE
        )
    if match is None:
        raise RuntimeError(
            f"Vina completed but no binding affinity could be parsed from {output_pdbqt}."
        )
    return float(match.group(1))


def run_vina(
    vina_executable,
    receptor_pdbqt,
    ligand_pdbqt,
    output_pdbqt,
    box,
    exhaustiveness=16,
    num_modes=10,
    seed=42,
):
    """
    Run AutoDock Vina using the supplied receptor, ligand and docking box.
    """

    print("\n[12] Running AutoDock Vina...")

    command = [
        vina_executable,

        "--receptor",
        receptor_pdbqt,

        "--ligand",
        ligand_pdbqt,

        "--center_x",
        str(box["center_x"]),

        "--center_y",
        str(box["center_y"]),

        "--center_z",
        str(box["center_z"]),

        "--size_x",
        str(box["size_x"]),

        "--size_y",
        str(box["size_y"]),

        "--size_z",
        str(box["size_z"]),

        "--exhaustiveness",
        str(exhaustiveness),

        "--num_modes",
        str(num_modes),

        "--seed",
        str(seed),

        "--out",
        output_pdbqt
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    print(result.stdout)

    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(
            "AutoDock Vina failed. "
            f"Exit code {result.returncode}. {result.stderr.strip()}"
        )

    if not Path(output_pdbqt).exists():
        raise RuntimeError("AutoDock Vina reported success but did not create its output file.")

    affinity = extract_best_affinity(output_pdbqt, result.stdout)

    print("Docking completed successfully.")
    print(f"Best binding affinity: {affinity:.2f} kcal/mol")
    print("Docking poses saved as:")
    print(output_pdbqt)

    return {
        "output_pdbqt": str(output_pdbqt),
        "binding_affinity": affinity,
        "seed": seed,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }
