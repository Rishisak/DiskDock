from Bio.PDB import PDBParser, PDBIO, Select
import subprocess


STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
}


class CleanProteinSelect(Select):

    def accept_residue(self, residue):

        # Keep only standard amino acids
        if residue.resname.strip() not in STANDARD_AMINO_ACIDS:
            return 0

        # Remove residues with negative residue numbers
        if residue.id[1] < 0:
            return 0

        return 1


def clean_protein(input_pdb, output_pdb):

    print("\n[9] Cleaning protein structure...")
    print("Input :", input_pdb)
    print("Output:", output_pdb)

    # Explicitly convert WindowsPath objects to strings
    input_pdb = str(input_pdb)
    output_pdb = str(output_pdb)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    io = PDBIO()
    io.set_structure(structure)

    # Save using a string filename
    io.save(
        output_pdb,
        CleanProteinSelect()
    )

    print("Protein cleaned successfully.")
    print("Saved as:", output_pdb)

    return output_pdb


def prepare_receptor(input_pdb, output_prefix):

    print("\n[10] Preparing protein receptor with Meeko...")
    print("Input PDB     :", input_pdb)
    print("Output prefix :", output_prefix)

    # Explicitly convert everything to strings
    input_pdb = str(input_pdb)
    output_prefix = str(output_prefix)

    command = [
        "mk_prepare_receptor.exe",
        "--read_pdb",
        input_pdb,
        "-o",
        output_prefix,
        "-p",

        # Use alternate location A when multiple locations exist
        "--default_altloc",
        "A",

        # Automatically remove residues that Meeko cannot
        # match to its standard residue templates
        "--allow_bad_res"
    ]

    print("\n----------------------------------------")
    print("Meeko command:")
    print("----------------------------------------")
    print(" ".join(command))
    print("----------------------------------------")

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True
        )

    except Exception as e:

        print("\n========================================")
        print("MEEKO EXECUTION ERROR")
        print("========================================")
        print("Could not execute mk_prepare_receptor.exe")
        print("Exception type:", type(e).__name__)
        print("Exception message:", str(e))
        print("========================================")

        raise RuntimeError(
            "Could not execute Meeko receptor preparation."
        ) from e

    print("\n========================================")
    print("MEEKO STDOUT")
    print("========================================")

    if result.stdout:
        print(result.stdout)
    else:
        print("(No stdout produced.)")

    print("\n========================================")
    print("MEEKO STDERR")
    print("========================================")

    if result.stderr:
        print(result.stderr)
    else:
        print("(No stderr produced.)")

    print("\n========================================")
    print("MEEKO RETURN CODE")
    print("========================================")
    print(result.returncode)

    if result.returncode != 0:

        print("\n========================================")
        print("MEEKO RECEPTOR PREPARATION FAILED")
        print("========================================")
        print("Input PDB:", input_pdb)
        print("Output prefix:", output_prefix)
        print("Return code:", result.returncode)

        raise RuntimeError(
            "Meeko receptor preparation failed. "
            "See the complete Meeko output above."
        )

    print("\nReceptor preparation completed successfully.")

    receptor_file = output_prefix + ".pdbqt"

    print("Saved as:", receptor_file)

    return receptor_file