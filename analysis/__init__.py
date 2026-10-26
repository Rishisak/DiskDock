"""Analysis package for DistDockNet.
Implements Step 12 of the pipeline: Target–Compound Intersection, PPI & Functional Analysis.
"""

from analysis.ppi import (
    fetch_string_ppi,
    build_ppi_graph,
    calculate_network_centralities,
    normalize_target_to_gene,
)
from analysis.enrichment import fetch_functional_enrichment
from analysis.intersection import build_target_compound_network, rank_compounds_by_polypharmacology
from analysis.visualization import generate_network_html, get_enrichment_altair_chart

__all__ = [
    "fetch_string_ppi",
    "build_ppi_graph",
    "calculate_network_centralities",
    "normalize_target_to_gene",
    "fetch_functional_enrichment",
    "build_target_compound_network",
    "rank_compounds_by_polypharmacology",
    "generate_network_html",
    "get_enrichment_altair_chart",
]
