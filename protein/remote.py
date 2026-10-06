"""Download and cache structures used by the docking batch."""

from pathlib import Path
import re

import requests
from Bio.PDB import MMCIFParser, PDBIO, Select


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = PROJECT_ROOT / "protein" / "cache"
RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"


class ProteinResidueSelect(Select):
    """Keep polymer residues when converting mmCIF to legacy PDB."""

    def accept_residue(self, residue):
        # Ligand IDs can exceed PDB's three-character residue-name field and
        # would corrupt the legacy fixed-column layout. The docking pipeline
        # removes hetero residues during cleaning anyway.
        return residue.id[0] == " "


def _safe_identifier(identifier, label):
    """Return an uppercase, filename-safe structure identifier."""
    value = identifier.strip().upper()
    if not re.fullmatch(r"[A-Z0-9]+", value):
        raise ValueError(f"Invalid {label}: {identifier!r}")
    return value


def _write_pdb_response(response, output_file, source_label):
    """Validate a downloaded PDB response before committing it to the cache."""
    if response.status_code != 200:
        raise RuntimeError(
            f"{source_label} download failed for {output_file.stem}. "
            f"HTTP {response.status_code}"
        )

    text = response.text
    if "ATOM  " not in text and "HETATM" not in text:
        raise RuntimeError(
            f"{source_label} returned no PDB atom records for {output_file.stem}."
        )

    temporary_file = output_file.with_suffix(".download")
    temporary_file.write_text(text, encoding="utf-8")
    temporary_file.replace(output_file)


def _write_mmcif_as_pdb(response, output_file, pdb_id):
    """Convert an RCSB mmCIF fallback to PDB for the downstream PDB pipeline."""
    if response.status_code != 200:
        raise RuntimeError(
            f"RCSB mmCIF download failed for {pdb_id}. HTTP {response.status_code}"
        )

    cif_file = output_file.with_suffix(".cif.download")
    pdb_file = output_file.with_suffix(".download")
    try:
        cif_file.write_text(response.text, encoding="utf-8")
        structure = MMCIFParser(QUIET=True).get_structure(pdb_id, str(cif_file))
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_file), ProteinResidueSelect())
        if not pdb_file.exists() or pdb_file.stat().st_size == 0:
            raise RuntimeError(f"RCSB mmCIF conversion produced no PDB atoms for {pdb_id}.")
        pdb_file.replace(output_file)
    except Exception as error:
        raise RuntimeError(f"Could not convert RCSB mmCIF structure {pdb_id} to PDB: {error}") from error
    finally:
        if cif_file.exists():
            cif_file.unlink()
        if pdb_file.exists():
            pdb_file.unlink()


def fetch_rcsb_structure(pdb_id):
    """Fetch an RCSB PDB structure once and return its local cached path."""
    pdb_id = _safe_identifier(pdb_id, "PDB ID")
    if len(pdb_id) != 4:
        raise ValueError(f"A PDB ID must contain four characters: {pdb_id!r}")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    output_file = CACHE_DIR / f"{pdb_id}.pdb"
    if output_file.exists() and output_file.stat().st_size > 0:
        print(f"[RCSB] Using cached structure: {output_file}")
        return str(output_file)

    print(f"[RCSB] Downloading {pdb_id}...")
    try:
        response = requests.get(RCSB_PDB_URL.format(pdb_id=pdb_id), timeout=60)
    except requests.RequestException as error:
        raise RuntimeError(f"RCSB download failed for {pdb_id}: {error}") from error

    if response.status_code == 404:
        print(f"[RCSB] PDB format unavailable for {pdb_id}; trying mmCIF...")
        try:
            cif_response = requests.get(RCSB_CIF_URL.format(pdb_id=pdb_id), timeout=60)
        except requests.RequestException as error:
            raise RuntimeError(f"RCSB mmCIF download failed for {pdb_id}: {error}") from error
        _write_mmcif_as_pdb(cif_response, output_file, pdb_id)
    else:
        _write_pdb_response(response, output_file, "RCSB")
    print(f"[RCSB] Saved: {output_file}")
    return str(output_file)


def fetch_alphafold_structure(uniprot_id):
    """Fetch the current AlphaFold DB PDB model for a UniProt accession.

    The AlphaFold API supplies the model-specific PDB URL, avoiding a hard-coded
    AlphaFold model-version suffix.
    """
    uniprot_id = _safe_identifier(uniprot_id, "UniProt ID")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    output_file = CACHE_DIR / f"AF-{uniprot_id}-F1.pdb"
    if output_file.exists() and output_file.stat().st_size > 0:
        print(f"[AlphaFold] Using cached structure: {output_file}")
        return str(output_file)

    print(f"[AlphaFold] Looking up {uniprot_id}...")
    try:
        metadata_response = requests.get(
            ALPHAFOLD_API_URL.format(uniprot_id=uniprot_id), timeout=60
        )
        if metadata_response.status_code != 200:
            raise RuntimeError(
                f"AlphaFold lookup failed for {uniprot_id}. "
                f"HTTP {metadata_response.status_code}"
            )
        predictions = metadata_response.json()
    except requests.RequestException as error:
        raise RuntimeError(f"AlphaFold lookup failed for {uniprot_id}: {error}") from error
    except ValueError as error:
        raise RuntimeError(
            f"AlphaFold returned invalid metadata for {uniprot_id}."
        ) from error

    if not predictions or not isinstance(predictions, list):
        raise RuntimeError(f"No AlphaFold prediction is available for {uniprot_id}.")

    pdb_url = predictions[0].get("pdbUrl")
    if not pdb_url:
        raise RuntimeError(
            f"AlphaFold metadata for {uniprot_id} does not include a PDB download URL."
        )

    print(f"[AlphaFold] Downloading {uniprot_id}...")
    try:
        pdb_response = requests.get(pdb_url, timeout=60)
    except requests.RequestException as error:
        raise RuntimeError(f"AlphaFold download failed for {uniprot_id}: {error}") from error

    _write_pdb_response(pdb_response, output_file, "AlphaFold")
    print(f"[AlphaFold] Saved: {output_file}")
    return str(output_file)


def fetch_structure(source, structure_id):
    """Fetch a structure from the supported source and return a local PDB path."""
    source = source.lower().strip()
    if source == "rcsb":
        return fetch_rcsb_structure(structure_id)
    if source == "alphafold":
        return fetch_alphafold_structure(structure_id)
    raise ValueError(f"Unknown structure source: {source}")
