import sqlite3
from pathlib import Path
import requests


BASE_DIR = Path(__file__).resolve().parent
DB_FILE = BASE_DIR / "distdocknet.db"


def get_connection():
    return sqlite3.connect(DB_FILE)


def initialize_database():
    connection = get_connection()
    cursor = connection.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS proteins (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pdb_id TEXT UNIQUE NOT NULL,
            pdb_data TEXT NOT NULL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS docking_runs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            workflow TEXT NOT NULL,
            protein_name TEXT NOT NULL,
            pdb_id TEXT,
            compound_name TEXT,
            smiles TEXT,
            reference_ligand TEXT,
            docking_mode TEXT,
            chain_id TEXT,
            exhaustiveness INTEGER,
            num_modes INTEGER,
            binding_affinity REAL,
            rmsd REAL,
            matched_atoms INTEGER,
            status TEXT NOT NULL,
            error_message TEXT,
            output_path TEXT,
            protein_path TEXT,
            ligand_path TEXT,
            receptor_path TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)

    connection.commit()
    connection.close()


def save_protein(pdb_id, pdb_data):
    connection = get_connection()
    cursor = connection.cursor()
    cursor.execute("""
        INSERT OR REPLACE INTO proteins
        (pdb_id, pdb_data)
        VALUES (?, ?)
    """, (pdb_id.upper(), pdb_data))
    connection.commit()
    connection.close()


def get_protein(pdb_id):
    connection = get_connection()
    cursor = connection.cursor()
    cursor.execute("""
        SELECT pdb_data
        FROM proteins
        WHERE pdb_id = ?
    """, (pdb_id.upper(),))
    row = cursor.fetchone()
    connection.close()

    if row is None:
        return None
    return row[0]


def get_or_fetch_protein(pdb_id):
    pdb_id = pdb_id.strip().upper()
    if not pdb_id:
        raise ValueError("PDB ID cannot be empty.")

    print(f"\nChecking database for PDB: {pdb_id}")
    pdb_data = get_protein(pdb_id)

    if pdb_data is not None:
        print(f"✓ {pdb_id} found in database.")
        return pdb_data

    print(f"✗ {pdb_id} not found in database.")
    print("Fetching from RCSB PDB...")

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print("URL:", url)

    try:
        response = requests.get(url, timeout=30)
    except requests.RequestException as error:
        raise RuntimeError(f"Could not connect to RCSB PDB: {error}")

    if response.status_code != 200:
        raise ValueError(
            f"PDB structure '{pdb_id}' could not be found on RCSB PDB."
        )

    pdb_data = response.text
    if not pdb_data.strip():
        raise ValueError(f"RCSB returned empty data for {pdb_id}.")

    save_protein(pdb_id, pdb_data)
    print(f"✓ {pdb_id} downloaded from RCSB.")
    print(f"✓ {pdb_id} saved to SQLite database.")
    return pdb_data


def create_working_pdb_file(pdb_id, output_directory):
    pdb_id = pdb_id.strip().upper()
    if not pdb_id:
        raise ValueError("PDB ID cannot be empty.")

    pdb_data = get_or_fetch_protein(pdb_id)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    pdb_file = output_directory / f"{pdb_id}.pdb"
    pdb_file.write_text(pdb_data, encoding="utf-8")

    print()
    print("----------------------------------------")
    print("WORKING PDB FILE")
    print("----------------------------------------")
    print("PDB ID :", pdb_id)
    print("File   :", pdb_file)
    print("----------------------------------------")

    return pdb_file


def save_docking_run(
    workflow,
    protein_name,
    pdb_id=None,
    compound_name=None,
    smiles=None,
    reference_ligand=None,
    docking_mode=None,
    chain_id=None,
    exhaustiveness=None,
    num_modes=None,
    binding_affinity=None,
    rmsd=None,
    matched_atoms=None,
    status="SUCCESS",
    error_message=None,
    output_path=None,
    protein_path=None,
    ligand_path=None,
    receptor_path=None,
):
    """Permanently save one docking/redocking run in SQLite."""
    connection = get_connection()
    cursor = connection.cursor()

    cursor.execute("""
        INSERT INTO docking_runs (
            workflow,
            protein_name,
            pdb_id,
            compound_name,
            smiles,
            reference_ligand,
            docking_mode,
            chain_id,
            exhaustiveness,
            num_modes,
            binding_affinity,
            rmsd,
            matched_atoms,
            status,
            error_message,
            output_path,
            protein_path,
            ligand_path,
            receptor_path
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        workflow,
        protein_name,
        pdb_id,
        compound_name,
        smiles,
        reference_ligand,
        docking_mode,
        chain_id,
        exhaustiveness,
        num_modes,
        binding_affinity,
        rmsd,
        matched_atoms,
        status,
        error_message,
        output_path,
        protein_path,
        ligand_path,
        receptor_path,
    ))

    run_id = cursor.lastrowid
    connection.commit()
    connection.close()
    return run_id


def get_docking_history(limit=500):
    """Return recent docking runs, newest first."""
    connection = get_connection()
    cursor = connection.cursor()

    cursor.execute("""
        SELECT
            id,
            created_at,
            workflow,
            protein_name,
            pdb_id,
            compound_name,
            reference_ligand,
            docking_mode,
            chain_id,
            exhaustiveness,
            num_modes,
            binding_affinity,
            rmsd,
            matched_atoms,
            status,
            error_message,
            output_path,
            protein_path,
            ligand_path,
            receptor_path,
            smiles
        FROM docking_runs
        ORDER BY id DESC
        LIMIT ?
    """, (int(limit),))

    rows = cursor.fetchall()
    columns = [column[0] for column in cursor.description]
    connection.close()

    return rows, columns


def get_docking_run(run_id):
    connection = get_connection()
    cursor = connection.cursor()

    cursor.execute("""
        SELECT *
        FROM docking_runs
        WHERE id = ?
    """, (int(run_id),))

    row = cursor.fetchone()
    columns = [column[0] for column in cursor.description]
    connection.close()

    if row is None:
        return None

    return dict(zip(columns, row))


initialize_database()
