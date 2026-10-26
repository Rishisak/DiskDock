"""Protein–Protein Interaction (PPI) Network Analysis Module.
Part of Step 12: Target–Compound Intersection, PPI & Functional Analysis.
Queries STRING DB API, constructs graph networks, and computes centrality/hub scores.
"""

from typing import List, Dict, Any, Optional
import requests
import networkx as nx
import pandas as pd
import numpy as np


STRING_API_NETWORK = "https://string-db.org/api/json/network"

GENE_ALIAS_MAP = {
    "TNF-ALPHA": "TNF",
    "TNFALPHA": "TNF",
    "IL-6": "IL6",
    "IL-1BETA": "IL1B",
    "IL1BETA": "IL1B",
    "IKKBETA": "IKBKB",
    "HIF-1ALPHA": "HIF1A",
    "HIF1ALPHA": "HIF1A",
    "MD-2": "LY96",
    "MD2": "LY96",
    "MMP-7": "MMP7",
    "MMP-13": "MMP13",
    "MMP-19": "MMP19",
    "ADAMTS-5": "ADAMTS5",
    "ADAMTS5": "ADAMTS5",
}


def normalize_target_to_gene(name: str) -> str:
    """Normalize a protein name or PDB label into a standard HGNC gene symbol.

    Examples:
      'TNF-alpha' -> 'TNF'
      'IL-6' -> 'IL6'
      'IKKbeta' -> 'IKBKB'
      'MMP-13' -> 'MMP13'
    """
    cleaned = str(name).strip().upper()
    if cleaned in GENE_ALIAS_MAP:
        return GENE_ALIAS_MAP[cleaned]
    # Remove hyphens/underscores if not in map
    simplified = cleaned.replace("-", "").replace("_", "")
    return GENE_ALIAS_MAP.get(simplified, simplified)


def fetch_string_ppi(
    genes: List[str],
    species: int = 9606,
    required_score: int = 400,
    timeout: int = 15,
) -> List[Dict[str, Any]]:
    """Fetch protein-protein interaction edges for a list of gene symbols from STRING DB.

    Args:
        genes: List of gene symbols (e.g., ['INS', 'AKT1', 'PPARG']).
        species: NCBI taxonomy ID (default: 9606 for Homo sapiens).
        required_score: Minimum interaction confidence score (0-1000, 400=medium, 700=high).
        timeout: Request timeout in seconds.

    Returns:
        List of interaction edge dictionaries from STRING DB.
    """
    cleaned_genes = [g.strip().upper() for g in genes if g and g.strip()]
    if not cleaned_genes:
        return []

    # STRING expects carriage-return separated identifiers
    params = {
        "identifiers": "\r".join(cleaned_genes),
        "species": species,
        "required_score": required_score,
        "caller_identity": "distdocknet_pipeline",
    }

    try:
        response = requests.post(
            STRING_API_NETWORK,
            data=params,
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()
        if isinstance(data, list):
            return data
        return []
    except Exception as exc:
        print(f"[analysis.ppi] STRING API error: {exc}")
        return []


def build_ppi_graph(
    edges: List[Dict[str, Any]],
    user_genes: Optional[List[str]] = None,
) -> nx.Graph:
    """Build a NetworkX Graph from STRING interaction records.

    Args:
        edges: Edge records from fetch_string_ppi.
        user_genes: Optional list of user genes to include even if they have 0 interactions.

    Returns:
        networkx.Graph with node & edge attributes.
    """
    G = nx.Graph()

    # Add isolated nodes if provided
    if user_genes:
        for gene in user_genes:
            cleaned = gene.strip().upper()
            if cleaned:
                G.add_node(cleaned, label=cleaned)

    for edge in edges:
        u = edge.get("preferredName_A")
        v = edge.get("preferredName_B")
        score = float(edge.get("score", 0.4))
        if not u or not v:
            continue

        u = u.strip().upper()
        v = v.strip().upper()

        G.add_node(u, label=u)
        G.add_node(v, label=v)
        G.add_edge(
            u,
            v,
            weight=score,
            score=score,
            escore=float(edge.get("escore", 0.0)),
            dscore=float(edge.get("dscore", 0.0)),
            tscore=float(edge.get("tscore", 0.0)),
        )

    return G


def calculate_network_centralities(
    G: nx.Graph,
) -> pd.DataFrame:
    """Compute centrality metrics and detect master Hub Targets in the PPI graph.

    Centrality metrics computed:
      - Degree (number of direct interactors)
      - Degree Centrality (normalized degree)
      - Betweenness Centrality (bottleneck / bridge metric)
      - Closeness Centrality (proximity to all other nodes)
      - Hub Score (composite normalized metric combining degree & betweenness)

    Args:
        G: NetworkX Graph.

    Returns:
        pd.DataFrame sorted descending by Hub Score.
    """
    if len(G.nodes) == 0:
        return pd.DataFrame(
            columns=[
                "gene",
                "degree",
                "degree_centrality",
                "betweenness_centrality",
                "closeness_centrality",
                "hub_score",
                "is_hub",
            ]
        )

    degrees = dict(G.degree())
    deg_centrality = nx.degree_centrality(G)

    # Betweenness centrality with edge weight inverted as distance if appropriate, or unweighted topology
    try:
        bet_centrality = nx.betweenness_centrality(G, normalized=True)
    except Exception:
        bet_centrality = {n: 0.0 for n in G.nodes}

    try:
        close_centrality = nx.closeness_centrality(G)
    except Exception:
        close_centrality = {n: 0.0 for n in G.nodes}

    rows = []
    for node in G.nodes:
        d = degrees.get(node, 0)
        dc = deg_centrality.get(node, 0.0)
        bc = bet_centrality.get(node, 0.0)
        cc = close_centrality.get(node, 0.0)
        rows.append(
            {
                "gene": node,
                "degree": d,
                "degree_centrality": round(dc, 4),
                "betweenness_centrality": round(bc, 4),
                "closeness_centrality": round(cc, 4),
            }
        )

    df = pd.DataFrame(rows)

    # Compute composite Hub Score: 0.6 * scaled(degree) + 0.4 * scaled(betweenness)
    max_d = df["degree"].max() if df["degree"].max() > 0 else 1.0
    max_bc = df["betweenness_centrality"].max() if df["betweenness_centrality"].max() > 0 else 1.0

    df["hub_score"] = (
        0.6 * (df["degree"] / max_d) + 0.4 * (df["betweenness_centrality"] / max_bc)
    ).round(4)

    # Mark as hub if in top 25% or degree >= mean + std
    mean_deg = df["degree"].mean()
    std_deg = df["degree"].std() if len(df) > 1 else 0
    threshold = mean_deg + (0.5 * std_deg) if std_deg > 0 else mean_deg

    df["is_hub"] = (df["degree"] >= threshold) & (df["degree"] > 1)

    df = df.sort_values(by=["hub_score", "degree"], ascending=[False, False]).reset_index(drop=True)
    return df
