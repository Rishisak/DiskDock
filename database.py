import sqlite3
from pathlib import Path
import requests
# ---------------------------------------------------------
# Database location
# ---------------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent

DB_FILE = BASE_DIR / "distdocknet.db"


# ---------------------------------------------------------
# Connect to database
# ---------------------------------------------------------

def get_connection():
    """
    Create a connection to the SQLite database.
    """

    return sqlite3.connect(DB_FILE)


# ---------------------------------------------------------
# Initialize database
# ---------------------------------------------------------

def initialize_database():
    """
    Create the required database tables if they
    do not already exist.
    """

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

    connection.commit()

    connection.close()


# ---------------------------------------------------------
# Save protein
# ---------------------------------------------------------

def save_protein(pdb_id, pdb_data):
    """
    Save a PDB structure into the database.
    """

    connection = get_connection()

    cursor = connection.cursor()

    cursor.execute("""
        INSERT OR REPLACE INTO proteins
        (pdb_id, pdb_data)
        VALUES (?, ?)
    """, (
        pdb_id.upper(),
        pdb_data
    ))

    connection.commit()

    connection.close()


# ---------------------------------------------------------
# Get protein
# ---------------------------------------------------------

def get_protein(pdb_id):
    """
    Retrieve a protein structure from the database.

    Returns:
        PDB text if found
        None if not found
    """

    connection = get_connection()

    cursor = connection.cursor()

    cursor.execute("""
        SELECT pdb_data
        FROM proteins
        WHERE pdb_id = ?
    """, (
        pdb_id.upper(),
    ))

    row = cursor.fetchone()

    connection.close()

    if row is None:
        return None

    return row[0]
def get_or_fetch_protein(pdb_id):
    """
    Get a protein from the SQLite database.

    If the protein is already stored in the database,
    return it directly.

    If it is not stored, download it from RCSB PDB,
    save it into the database, and return it.
    """

    # Clean the PDB ID
    pdb_id = pdb_id.strip().upper()

    if not pdb_id:
        raise ValueError("PDB ID cannot be empty.")

    # -----------------------------------------------------
    # Step 1: Check SQLite database
    # -----------------------------------------------------

    print(f"\nChecking database for PDB: {pdb_id}")

    pdb_data = get_protein(pdb_id)

    if pdb_data is not None:

        print(f"✓ {pdb_id} found in database.")

        return pdb_data

    # -----------------------------------------------------
    # Step 2: Protein not found
    # -----------------------------------------------------

    print(f"✗ {pdb_id} not found in database.")

    print("Fetching from RCSB PDB...")

    # -----------------------------------------------------
    # Step 3: RCSB download URL
    # -----------------------------------------------------

    url = (
        f"https://files.rcsb.org/download/"
        f"{pdb_id}.pdb"
    )

    print("URL:", url)

    # -----------------------------------------------------
    # Step 4: Download
    # -----------------------------------------------------

    try:

        response = requests.get(
            url,
            timeout=30
        )

    except requests.RequestException as error:

        raise RuntimeError(
            f"Could not connect to RCSB PDB: {error}"
        )

    # -----------------------------------------------------
    # Step 5: Check response
    # -----------------------------------------------------

    if response.status_code != 200:

        raise ValueError(
            f"PDB structure '{pdb_id}' "
            f"could not be found on RCSB PDB."
        )

    pdb_data = response.text

    # -----------------------------------------------------
    # Step 6: Make sure something was downloaded
    # -----------------------------------------------------

    if not pdb_data.strip():

        raise ValueError(
            f"RCSB returned empty data for {pdb_id}."
        )

    # -----------------------------------------------------
    # Step 7: Save protein in SQLite
    # -----------------------------------------------------

    save_protein(
        pdb_id,
        pdb_data
    )

    print(f"✓ {pdb_id} downloaded from RCSB.")

    print(f"✓ {pdb_id} saved to SQLite database.")

    # -----------------------------------------------------
    # Step 8: Return PDB data
    # -----------------------------------------------------

    return pdb_data
def create_working_pdb_file(pdb_id, output_directory):
    """
    Get a protein from the SQLite database/RCSB
    and create a temporary PDB file for docking.
    """

    pdb_id = pdb_id.strip().upper()

    if not pdb_id:
        raise ValueError("PDB ID cannot be empty.")

    # Get PDB data from database or RCSB
    pdb_data = get_or_fetch_protein(pdb_id)

    # Create output directory
    output_directory = Path(output_directory)

    output_directory.mkdir(
        parents=True,
        exist_ok=True
    )

    # Create temporary working PDB file
    pdb_file = output_directory / f"{pdb_id}.pdb"

    # Write PDB data
    pdb_file.write_text(
        pdb_data,
        encoding="utf-8"
    )

    print()
    print("----------------------------------------")
    print("WORKING PDB FILE")
    print("----------------------------------------")
    print("PDB ID :", pdb_id)
    print("File   :", pdb_file)
    print("----------------------------------------")

    return pdb_file