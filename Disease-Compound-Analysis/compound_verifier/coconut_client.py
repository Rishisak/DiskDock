"""
coconut_client.py
-----------------
COCONUT (COlleCtion of Open Natural prodUcTs) REST API client for Step 7 verification.

Validates chemical natural origin against 400,000+ curated natural products via
their Standard InChIKey, retrieving authoritative COCONUT identifiers (CNPxxxxx),
natural source collections, and organism counts.

Integrated with SQLite persistent caching and graceful failover.
COCONUT 2.0 Web Portal: https://coconut.naturalproducts.net
"""

import logging
import re
from typing import Any, Dict, Optional
import requests

try:
    from compound_cache import get_generic_cache, set_generic_cache
except ImportError:
    try:
        from ..compound_cache import get_generic_cache, set_generic_cache
    except Exception:
        def get_generic_cache(ns, k): return None
        def set_generic_cache(ns, k, v): pass

logger = logging.getLogger(__name__)

COCONUT_SEARCH_URL = "https://coconut.naturalproducts.net/api/search"


class CoconutClient:
    """Client for the open COCONUT 2.0 natural products database with persistent caching."""

    def __init__(self, timeout: float = 3.5):
        self.session = requests.Session()
        self.session.headers.update({
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "Biomedical-AI-Compound-Verifier/1.0 (academic/research)",
        })
        self.timeout = timeout

    def lookup_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Query COCONUT for an exact Standard InChIKey match.
        Returns parsed natural origin evidence metadata or None if not found/offline.
        """
        clean_key = (inchikey or "").strip().upper()
        if not clean_key or not self.is_valid_inchikey(clean_key):
            return None

        # Check persistent SQLite cache first (0ms latency)
        cached = get_generic_cache("coconut", clean_key)
        if cached is not None:
            return cached if cached != "__NOT_FOUND__" else None

        payload = {"query": clean_key}

        try:
            resp = self.session.post(
                COCONUT_SEARCH_URL,
                json=payload,
                timeout=self.timeout,
            )
            if resp.status_code == 200:
                res_json = resp.json()
                items = res_json.get("data", {}).get("data", [])
                if items:
                    mol = items[0]
                    result = {
                        "coconut_id": mol.get("identifier"),
                        "name": mol.get("name"),
                        "canonical_smiles": mol.get("canonical_smiles"),
                        "collection_count": mol.get("collection_count", 0),
                        "organism_count": mol.get("organism_count", 0),
                        "citation_count": mol.get("citation_count", 0),
                        "annotation_level": mol.get("annotation_level", 0),
                        "is_natural": True,
                    }
                    set_generic_cache("coconut", clean_key, result)
                    return result

                # 200 OK but empty result
                set_generic_cache("coconut", clean_key, "__NOT_FOUND__")
                return None
            else:
                logger.debug(f"COCONUT returned HTTP {resp.status_code} for {clean_key}")
                return None
        except Exception as exc:
            logger.debug(f"COCONUT lookup skipped/timed out for {clean_key}: {exc}")
            return None

    @staticmethod
    def is_valid_inchikey(key: Optional[str]) -> bool:
        """Validate whether a string conforms to Standard InChIKey syntax (14-10-1 format)."""
        if not key or not isinstance(key, str):
            return False
        pattern = r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$"
        return bool(re.match(pattern, key.strip()))

