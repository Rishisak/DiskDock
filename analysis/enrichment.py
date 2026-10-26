"""Functional and Pathway Enrichment Module.
Part of Step 12: Target–Compound Intersection, PPI & Functional Analysis.
Queries STRING Enrichment API for KEGG pathways, GO terms (BP, MF, CC), and Reactome pathways.
"""

from typing import List, Dict, Any, Optional
import math
import requests
import pandas as pd


STRING_API_ENRICHMENT = "https://string-db.org/api/json/enrichment"

CATEGORY_MAPPING = {
    "KEGG": "KEGG Pathways",
    "Process": "GO Biological Process",
    "Function": "GO Molecular Function",
    "COMPARTMENTS": "GO Cellular Component",
    "RCTM": "Reactome Pathways",
    "WikiPathways": "WikiPathways",
    "DISEASES": "Disease Associations",
}


def fetch_functional_enrichment(
    genes: List[str],
    species: int = 9606,
    fdr_threshold: float = 0.05,
    timeout: int = 15,
) -> Dict[str, pd.DataFrame]:
    """Fetch functional enrichment results for a gene set across KEGG, GO, and Reactome.

    Args:
        genes: List of gene symbols (e.g., ['INS', 'AKT1', 'PPARG']).
        species: NCBI taxonomy ID (9606 = Homo sapiens).
        fdr_threshold: Maximum false discovery rate (FDR) cutoff.
        timeout: Request timeout in seconds.

    Returns:
        Dictionary mapping category names (e.g., 'KEGG Pathways', 'GO Biological Process')
        to DataFrames sorted by statistical significance (lowest FDR).
    """
    cleaned_genes = [g.strip().upper() for g in genes if g and g.strip()]
    empty_result = {v: pd.DataFrame() for v in CATEGORY_MAPPING.values()}

    if not cleaned_genes:
        return empty_result

    params = {
        "identifiers": "\r".join(cleaned_genes),
        "species": species,
        "caller_identity": "distdocknet_pipeline",
    }

    try:
        response = requests.post(
            STRING_API_ENRICHMENT,
            data=params,
            timeout=timeout,
        )
        response.raise_for_status()
        records = response.json()
    except Exception as exc:
        print(f"[analysis.enrichment] STRING API error: {exc}")
        return empty_result

    if not isinstance(records, list):
        return empty_result

    categorized_data: Dict[str, List[Dict[str, Any]]] = {
        label: [] for label in CATEGORY_MAPPING.values()
    }

    for item in records:
        raw_cat = item.get("category", "")
        if raw_cat not in CATEGORY_MAPPING:
            continue

        label = CATEGORY_MAPPING[raw_cat]
        fdr = float(item.get("fdr", 1.0))
        if fdr > fdr_threshold:
            continue

        # Compute -log10(FDR)
        neg_log10_fdr = -math.log10(max(fdr, 1e-30))

        input_genes = item.get("inputGenes", [])
        if isinstance(input_genes, list):
            gene_str = ", ".join(input_genes)
            gene_count = len(input_genes)
        else:
            gene_str = str(input_genes)
            gene_count = int(item.get("number_of_genes", 0))

        categorized_data[label].append(
            {
                "term": item.get("term", ""),
                "description": item.get("description", ""),
                "fdr": fdr,
                "neg_log10_fdr": round(neg_log10_fdr, 2),
                "gene_count": gene_count,
                "genes": gene_str,
            }
        )

    dfs: Dict[str, pd.DataFrame] = {}
    for label, items in categorized_data.items():
        if items:
            df = pd.DataFrame(items)
            df = df.sort_values(by=["fdr", "gene_count"], ascending=[True, False]).reset_index(drop=True)
            dfs[label] = df
        else:
            dfs[label] = pd.DataFrame(
                columns=["term", "description", "fdr", "neg_log10_fdr", "gene_count", "genes"]
            )

    return dfs
