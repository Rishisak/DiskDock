"""
verifier.py
-----------
Master Orchestrator for Step 7: Compound Identity & Natural-Origin Verification.

Coordinates:
  1. Chemical Identity Resolution via PubChem PUG REST API
     - Resolves canonical name, IUPAC name, and chemical synonyms
     - Standardizes SMILES and Standard InChIKey (14-10-1 format)
     - Retrieves authoritative PubChem CID
  2. Natural vs. Synthetic Origin Classification
     - Cross-references curated phytochemical registries and synthetic drugs
     - Cross-references COCONUT 2.0 open repository (400,000+ natural products)
     - Incorporates ChEMBL 'natural_product' flags and PubChem xrefs
     - Assigns confidence scores [0.0 - 1.0] and phytochemical taxonomy
  3. Provenance Retention
     - Preserves all Step 6 retrieval routes, target relations, and bioactivities
     - Appends full audit trail of chemical resolution and verification metrics
  4. Concurrent Batch Processing
     - Evaluates multiple candidates concurrently with thread pool execution
     - Integrates with persistent SQLite caching for zero-delay repeated lookups
"""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import logging
from typing import Any, Dict, List, Optional, Union

from .coconut_client import CoconutClient
from .models import (
    CompoundIdentity,
    CompoundVerificationBatchResult,
    NaturalOriginStatus,
    NaturalOriginVerification,
    VerifiedCompound,
)
from .natural_classifier import NaturalOriginClassifier
from .pubchem_client import PubChemClient

logger = logging.getLogger(__name__)


class CompoundVerifier:
    """
    Step 7 Engine: Standardizes compound chemical identities, verifies
    natural origins, distinguishes natural vs. synthetic compounds, and
    retains comprehensive evidence provenance.
    """

    def __init__(self, timeout: int = 15, enable_pubchem: bool = True, enable_coconut: bool = True):
        self.enable_pubchem = enable_pubchem
        self.pubchem = PubChemClient(timeout=timeout) if enable_pubchem else None
        self.enable_coconut = enable_coconut
        self.coconut = CoconutClient(timeout=3.5) if enable_coconut else None
        self.classifier = NaturalOriginClassifier()

    def verify_compound(
        self,
        candidate: Union[Dict[str, Any], Any],
    ) -> VerifiedCompound:
        """
        Verify an individual compound candidate (CompoundCandidate model or dict).
        Resolves identity, determines natural vs synthetic status, and retains provenance.
        """
        # 1. Normalize input candidate attributes
        if isinstance(candidate, dict):
            c_name = candidate.get("name", "").strip()
            c_id = candidate.get("compound_id") or candidate.get("chembl_id") or f"CMPD:{c_name.upper()}"
            c_synonyms = list(candidate.get("synonyms") or [])
            c_smiles = candidate.get("smiles")
            c_inchikey = candidate.get("inchikey")
            c_pubchem_cid = candidate.get("pubchem_cid")
            c_chembl_id = candidate.get("chembl_id")
            c_routes = list(candidate.get("retrieval_routes") or [])
            c_targets = list(candidate.get("targets_targeted") or [])
            c_target_relations = list(candidate.get("target_relations") or [])
            c_sources = list(candidate.get("evidence_sources") or [])
            c_citations = int(candidate.get("citation_count") or 0)
            c_potency = candidate.get("best_bioactivity_nM")
            c_provenance = dict(candidate.get("provenance") or {})
            c_natural_hint = candidate.get("natural_origin")
            c_chembl_meta = c_provenance.get("chembl_meta")
        else:
            # Pydantic CompoundCandidate instance
            c_name = getattr(candidate, "name", "").strip()
            c_id = getattr(candidate, "compound_id", "") or f"CMPD:{c_name.upper()}"
            c_synonyms = list(getattr(candidate, "synonyms", []))
            c_smiles = getattr(candidate, "smiles", None)
            c_inchikey = getattr(candidate, "inchikey", None)
            c_pubchem_cid = getattr(candidate, "pubchem_cid", None)
            c_chembl_id = getattr(candidate, "chembl_id", None)
            c_routes = list(getattr(candidate, "retrieval_routes", []))
            c_targets = list(getattr(candidate, "targets_targeted", []))
            c_target_relations = list(getattr(candidate, "target_relations", []))
            c_sources = list(getattr(candidate, "evidence_sources", []))
            c_citations = int(getattr(candidate, "citation_count", 0))
            c_potency = getattr(candidate, "best_bioactivity_nM", None)
            c_provenance = dict(getattr(candidate, "provenance", {}))
            c_natural_hint = getattr(candidate, "natural_origin", None)
            c_chembl_meta = c_provenance.get("chembl_meta")

        # 2. Resolve chemical identity via PubChem
        resolved_identity_dict: Optional[Dict[str, Any]] = None
        if self.pubchem and self.enable_pubchem:
            try:
                resolved_identity_dict = self.pubchem.resolve_identity(
                    name=c_name,
                    inchikey=c_inchikey,
                    smiles=c_smiles,
                    pubchem_cid=c_pubchem_cid,
                )
            except Exception as exc:
                logger.warning(f"Error during PubChem identity resolution for '{c_name}': {exc}")

        if resolved_identity_dict:
            if not resolved_identity_dict.get("name"):
                resolved_identity_dict["name"] = c_name
            identity = CompoundIdentity(**resolved_identity_dict)
            identity_resolved = True
        else:
            # Fallback to local data provided in candidate
            identity = CompoundIdentity(
                name=c_name,
                pubchem_cid=c_pubchem_cid,
                inchikey=c_inchikey,
                canonical_smiles=c_smiles,
                isomeric_smiles=None,
                iupac_name=None,
                molecular_formula=None,
                molecular_weight=None,
                synonyms=c_synonyms[:10],
            )
            identity_resolved = False

        # 3. Check COCONUT natural products database via InChIKey
        coconut_data = None
        if identity.inchikey and self.coconut and self.enable_coconut:
            try:
                coconut_data = self.coconut.lookup_inchikey(identity.inchikey)
            except Exception as exc:
                logger.debug(f"COCONUT lookup exception for {c_name}: {exc}")

        # 4. Verify natural origin & distinguish natural vs synthetic
        natural_verification = self.classifier.classify(
            name=c_name,
            synonyms=identity.synonyms,
            chembl_meta=c_chembl_meta,
            retrieval_routes=c_routes,
            is_chembl_natural=c_natural_hint,
            coconut_data=coconut_data,
        )

        # 5. Determine final persistent compound ID
        # Prefer PubChem CID if resolved, else ChEMBL ID, else original candidate ID
        if identity.pubchem_cid:
            final_compound_id = f"CID:{identity.pubchem_cid}"
        elif c_chembl_id:
            final_compound_id = c_chembl_id
        else:
            final_compound_id = c_id

        # 6. Build enriched provenance retaining all Step 6 data
        verification_provenance = {
            "step_7_verified_at": datetime.now(timezone.utc).isoformat(),
            "pubchem_resolution": "success" if identity_resolved else "unmatched_or_offline",
            "pubchem_cid": identity.pubchem_cid,
            "inchikey_standardized": bool(identity.inchikey and PubChemClient.is_valid_inchikey(identity.inchikey)),
            "origin_status": natural_verification.status.value,
            "origin_confidence": natural_verification.confidence,
            "determination_method": natural_verification.determination_method,
            "evidence_sources": natural_verification.evidence_sources,
            "coconut_id": natural_verification.coconut_id,
            "coconut_collections": natural_verification.coconut_collections,
            "coconut_organisms": natural_verification.coconut_organisms,
        }
        # Merge with existing Step 6 provenance without overwriting original PMIDs or snippets
        merged_provenance = {**c_provenance, "verification": verification_provenance}

        # Update evidence sources list
        extra_sources = []
        if identity_resolved:
            extra_sources.append("PubChem")
        if coconut_data:
            extra_sources.append("COCONUT")
        updated_sources = list(set(c_sources + extra_sources))

        # Determine dual-route hit flag
        is_dual = (
            "disease_natural_literature" in c_routes
            and ("protein_chembl_bioactivity" in c_routes or "protein_literature" in c_routes)
        )

        return VerifiedCompound(
            compound_id=final_compound_id,
            name=c_name.lower(),
            identity=identity,
            natural_verification=natural_verification,
            retrieval_routes=c_routes,
            targets_targeted=c_targets,
            target_relations=c_target_relations,
            evidence_sources=updated_sources,
            citation_count=c_citations,
            best_bioactivity_nM=c_potency,
            is_dual_route=is_dual,
            provenance=merged_provenance,
        )

    def verify_batch(
        self,
        candidates: List[Union[Dict[str, Any], Any]],
        disease_name: Optional[str] = None,
    ) -> CompoundVerificationBatchResult:
        """
        Verify a batch of candidate compounds from Step 6 concurrently.
        Deduplicates, standardizes, classifies, and sorts results by natural product validity.
        """
        verified_compounds: List[VerifiedCompound] = []
        seen_keys = set()

        # Concurrent resolution with thread pool
        if len(candidates) > 1:
            max_workers = min(8, len(candidates))
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                raw_verified = list(executor.map(self.verify_compound, candidates))
        else:
            raw_verified = [self.verify_compound(c) for c in candidates]

        for verified in raw_verified:
            # Deduplication key preference: InChIKey > PubChem CID > Normalized Name
            if verified.identity.inchikey:
                dedup_key = f"IK:{verified.identity.inchikey}"
            elif verified.identity.pubchem_cid:
                dedup_key = f"CID:{verified.identity.pubchem_cid}"
            else:
                dedup_key = f"NAME:{verified.name.strip().lower()}"

            if dedup_key not in seen_keys:
                seen_keys.add(dedup_key)
                verified_compounds.append(verified)
            else:
                # Merge target relations and citations into existing entry
                for existing in verified_compounds:
                    e_key = (
                        f"IK:{existing.identity.inchikey}"
                        if existing.identity.inchikey
                        else f"CID:{existing.identity.pubchem_cid}"
                        if existing.identity.pubchem_cid
                        else f"NAME:{existing.name.strip().lower()}"
                    )
                    if e_key == dedup_key:
                        # Combine target symbols
                        for sym in verified.targets_targeted:
                            if sym not in existing.targets_targeted:
                                existing.targets_targeted.append(sym)
                        # Combine retrieval routes
                        for route in verified.retrieval_routes:
                            if route not in existing.retrieval_routes:
                                existing.retrieval_routes.append(route)
                        # Re-evaluate dual route flag
                        existing.is_dual_route = (
                            "disease_natural_literature" in existing.retrieval_routes
                            and (
                                "protein_chembl_bioactivity" in existing.retrieval_routes
                                or "protein_literature" in existing.retrieval_routes
                            )
                        )
                        # Update best bioactivity
                        if verified.best_bioactivity_nM is not None:
                            if (
                                existing.best_bioactivity_nM is None
                                or verified.best_bioactivity_nM < existing.best_bioactivity_nM
                            ):
                                existing.best_bioactivity_nM = verified.best_bioactivity_nM
                        existing.citation_count += verified.citation_count
                        break

        # Calculate breakdown counts
        natural_count = sum(
            1 for c in verified_compounds
            if c.natural_verification.status == NaturalOriginStatus.NATURAL
        )
        derivative_count = sum(
            1 for c in verified_compounds
            if c.natural_verification.status == NaturalOriginStatus.NATURAL_DERIVATIVE
        )
        synthetic_count = sum(
            1 for c in verified_compounds
            if c.natural_verification.status == NaturalOriginStatus.SYNTHETIC
        )
        unverified_count = sum(
            1 for c in verified_compounds
            if c.natural_verification.status == NaturalOriginStatus.UNVERIFIED
        )

        # Sort:
        # 1. Natural / Natural Derivative first
        # 2. Dual-route priority
        # 3. Known bioactivity
        # 4. Citations count (descending)
        def _sort_key(c: VerifiedCompound):
            status_order = {
                NaturalOriginStatus.NATURAL: 0,
                NaturalOriginStatus.NATURAL_DERIVATIVE: 1,
                NaturalOriginStatus.UNVERIFIED: 2,
                NaturalOriginStatus.SYNTHETIC: 3,
            }
            s_rank = status_order.get(c.natural_verification.status, 2)
            d_rank = 0 if c.is_dual_route else 1
            has_bio = 0 if c.best_bioactivity_nM is not None else 1
            return (
                s_rank,
                d_rank,
                has_bio,
                -c.citation_count,
                c.name,
            )

        sorted_compounds = sorted(verified_compounds, key=_sort_key)

        return CompoundVerificationBatchResult(
            disease_name=disease_name,
            total_candidates=len(candidates),
            verified_natural_count=natural_count,
            natural_derivative_count=derivative_count,
            synthetic_count=synthetic_count,
            unverified_count=unverified_count,
            compounds=sorted_compounds,
        )
