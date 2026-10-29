"""
disease_docking_bridge.py - Bridge connecting Disease-Compound-Analysis
with the DistDockNet Streamlit Molecular Docking & Step 12 Systems Pharmacology Pipeline.
"""

import json
import logging
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import requests

logger = logging.getLogger(__name__)

BACKEND_URL = "http://localhost:8000"

# Curated High-Affinity / Co-Crystallized Structure PDB Mapping
# Provides instant, zero-latency, high-resolution structures for top therapeutic targets
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
}


def is_backend_healthy(timeout: float = 1.5) -> bool:
    """Check if the Disease-Compound-Analysis server is listening on port 8000."""
    try:
        resp = requests.get(f"{BACKEND_URL}/api/health", timeout=timeout)
        return resp.status_code == 200
    except Exception:
        return False


def get_disease_suggestions(query: str, limit: int = 8) -> List[Dict[str, Any]]:
    """Fetch disease name autocomplete suggestions from the backend normalizer."""
    if not query or len(query.strip()) < 2:
        return []
    try:
        url = f"{BACKEND_URL}/api/suggest?q={urllib.parse.quote(query.strip())}"
        resp = requests.get(url, timeout=3.0)
        if resp.status_code == 200:
            data = resp.json()
            return data.get("suggestions", [])[:limit]
    except Exception as exc:
        logger.warning(f"Failed to fetch disease suggestions: {exc}")
    return []


def run_disease_target_discovery(
    disease_query: str,
    target_limit: int = 6,
    papers_per_target: int = 2,
    min_score: float = 0.40,
    timeout: float = 45.0,
) -> Dict[str, Any]:
    """
    Execute disease normalization, Open Targets querying, and PubMed evidence fusion.
    Returns normalized disease metadata, synonyms, and ranked targets.
    """
    payload = {
        "disease": disease_query.strip(),
        "target_limit": int(target_limit),
        "papers_per_target": int(papers_per_target),
        "min_score": float(min_score),
    }
    resp = requests.post(
        f"{BACKEND_URL}/api/analyze",
        json=payload,
        timeout=timeout,
    )
    if resp.status_code != 200:
        err_msg = resp.text
        try:
            err_json = resp.json()
            err_msg = err_json.get("message", err_msg)
        except Exception:
            pass
        raise RuntimeError(f"Target analysis failed: {err_msg}")

    data = resp.json()
    disease_obj = data.get("disease", {}) if isinstance(data.get("disease"), dict) else {}
    canonical_name = disease_obj.get("canonical_name") or data.get("disease_name") or disease_query

    return {
        "status": data.get("status", "success"),
        "disease_name": canonical_name,
        "disease_id": disease_obj.get("canonical_id", ""),
        "synonyms": disease_obj.get("synonyms", []),
        "targets": data.get("targets", []),
        "raw": data,
    }


def run_natural_compound_discovery(
    disease_name: str,
    targets: List[Dict[str, Any]],
    synonyms: Optional[List[str]] = None,
    timeout: float = 60.0,
) -> Dict[str, Any]:
    """
    Execute Route A (literature) + Route B (ChEMBL bioactivities),
    followed by PubChem identity resolution and COCONUT database verification.
    """
    # Clean targets to pass to /api/compounds
    clean_targets = []
    for t in targets:
        clean_targets.append({
            "target_symbol": t.get("target_symbol", ""),
            "target_name": t.get("target_name", ""),
            "fused_score": t.get("fused_score", 0.0),
        })

    payload = {
        "disease": disease_name.strip(),
        "targets": clean_targets,
        "synonyms": synonyms or [],
    }
    resp = requests.post(
        f"{BACKEND_URL}/api/compounds",
        json=payload,
        timeout=timeout,
    )
    if resp.status_code != 200:
        err_msg = resp.text
        try:
            err_json = resp.json()
            err_msg = err_json.get("message", err_msg)
        except Exception:
            pass
        raise RuntimeError(f"Compound discovery failed: {err_msg}")

    return resp.json()


def resolve_target_pdb(target_symbol: str, target_name: str = "") -> Tuple[str, str]:
    """
    Resolve a gene symbol / target protein to a 4-letter experimental PDB ID.
    Returns (pdb_id, source_method).
    """
    sym = (target_symbol or "").strip().upper()
    if sym in CURATED_TARGET_PDBS:
        return CURATED_TARGET_PDBS[sym], "Curated High-Res Crystal Structure"

    # Query RCSB PDB using search_pdb from protein.lookup
    try:
        from protein.lookup import search_pdb
        pdb_id = search_pdb(sym)
        if pdb_id and len(pdb_id) == 4:
            return pdb_id.upper(), "RCSB PDB Dynamic Query"
    except Exception as exc:
        logger.warning(f"RCSB search for symbol {sym} failed: {exc}")

    # Fallback to searching with target name if available
    if target_name and target_name.strip():
        try:
            from protein.lookup import search_pdb
            pdb_id = search_pdb(target_name.strip())
            if pdb_id and len(pdb_id) == 4:
                return pdb_id.upper(), "RCSB PDB Dynamic Query (Name)"
        except Exception:
            pass

    return "", "Unresolved"
