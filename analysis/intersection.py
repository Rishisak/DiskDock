"""Target–Compound Intersection & Polypharmacology Analysis Module.
Part of Step 12: Target–Compound Intersection, PPI & Functional Analysis.
Analyzes binding affinity matrices, constructs bipartite interaction networks,
and ranks multi-target nutraceutical compounds targeting disease hub proteins.
"""

from typing import List, Dict, Any, Optional
import pandas as pd
import numpy as np
import networkx as nx

from analysis.ppi import normalize_target_to_gene


def build_target_compound_network(
    docking_df: pd.DataFrame,
    affinity_threshold: float = -7.0,
    target_col: str = "protein_name",
    compound_col: str = "compound_name",
    affinity_col: str = "binding_affinity",
) -> nx.Graph:
    """Construct a bipartite NetworkX graph connecting compounds to their targeted proteins.

    Args:
        docking_df: DataFrame with docking results.
        affinity_threshold: Cutoff (kcal/mol); interactions with affinity <= threshold are included.
        target_col: Column name for protein/target identifier.
        compound_col: Column name for compound identifier.
        affinity_col: Column name for binding affinity in kcal/mol.

    Returns:
        networkx.Graph containing bipartite nodes and weighted edges.
    """
    G = nx.Graph()

    if docking_df.empty:
        return G

    # Standardize column lookups
    cols = {c.lower(): c for c in docking_df.columns}
    t_col = cols.get(target_col.lower(), target_col)
    c_col = cols.get(compound_col.lower(), compound_col)
    a_col = cols.get(affinity_col.lower(), affinity_col)

    if t_col not in docking_df or c_col not in docking_df or a_col not in docking_df:
        return G

    # Filter by affinity threshold (more negative = stronger binding)
    valid_hits = docking_df[docking_df[a_col] <= affinity_threshold].copy()

    for _, row in valid_hits.iterrows():
        p_name = normalize_target_to_gene(row[t_col])
        c_name = str(row[c_col]).strip()
        aff = float(row[a_col])

        if not p_name or not c_name or np.isnan(aff):
            continue

        # Add target node
        G.add_node(
            p_name,
            node_type="target",
            label=p_name,
            bipartite=0,
        )

        # Add compound node
        G.add_node(
            c_name,
            node_type="compound",
            label=c_name,
            bipartite=1,
        )

        # Edge weight proportional to binding strength: weight = abs(affinity)
        G.add_edge(
            c_name,
            p_name,
            affinity=aff,
            weight=abs(aff),
        )

    return G


def rank_compounds_by_polypharmacology(
    docking_df: pd.DataFrame,
    hub_df: Optional[pd.DataFrame] = None,
    affinity_threshold: float = -7.0,
    target_col: str = "protein_name",
    compound_col: str = "compound_name",
    affinity_col: str = "binding_affinity",
) -> pd.DataFrame:
    """Rank compounds by multi-target polypharmacology and affinity against disease hub proteins.

    Compounds that bind multiple disease-relevant targets with high affinity are prioritized.

    Args:
        docking_df: Docking run results DataFrame.
        hub_df: Optional hub targets DataFrame from calculate_network_centralities.
        affinity_threshold: Cutoff in kcal/mol for strong binding.
        target_col: Column name for target protein.
        compound_col: Column name for compound name.
        affinity_col: Column name for binding affinity.

    Returns:
        pd.DataFrame ranked by polypharmacology score.
    """
    if docking_df.empty:
        return pd.DataFrame()

    cols = {c.lower(): c for c in docking_df.columns}
    t_col = cols.get(target_col.lower(), target_col)
    c_col = cols.get(compound_col.lower(), compound_col)
    a_col = cols.get(affinity_col.lower(), affinity_col)

    if t_col not in docking_df or c_col not in docking_df or a_col not in docking_df:
        return pd.DataFrame()

    df = docking_df.copy()
    df[t_col] = df[t_col].apply(normalize_target_to_gene)
    df[c_col] = df[c_col].astype(str).str.strip()
    df[a_col] = pd.to_numeric(df[a_col], errors="coerce")
    df = df.dropna(subset=[a_col])

    # Extract hub genes set if provided
    hub_genes = set()
    if hub_df is not None and not hub_df.empty and "gene" in hub_df.columns:
        hub_genes = set(hub_df[hub_df.get("is_hub", True)]["gene"].astype(str).str.upper())

    results = []
    for compound, group in df.groupby(c_col):
        total_targets = group[t_col].nunique()
        strong_hits = group[group[a_col] <= affinity_threshold]
        strong_count = strong_hits[t_col].nunique()

        best_affinity = group[a_col].min()
        mean_affinity = round(group[a_col].mean(), 2)

        # Targeted proteins list
        strong_targets = sorted(list(strong_hits[t_col].unique()))

        # Check hub targets hit
        hit_hubs = [t for t in strong_targets if t in hub_genes]
        hub_count = len(hit_hubs)

        # Polypharmacology Score calculation:
        # Rewards:
        #   - Number of strong target hits (strong_count * 2.0)
        #   - Strength of best affinity (abs(best_affinity) * 0.5)
        #   - Bonus for hitting master disease hubs (hub_count * 3.0)
        poly_score = (
            (strong_count * 2.0)
            + (abs(min(best_affinity, 0.0)) * 0.5)
            + (hub_count * 3.0)
        )

        results.append(
            {
                "compound": compound,
                "poly_score": round(poly_score, 2),
                "strong_targets_count": strong_count,
                "total_targets_tested": total_targets,
                "best_affinity_kcal_mol": round(best_affinity, 2),
                "mean_affinity_kcal_mol": mean_affinity,
                "targeted_proteins": ", ".join(strong_targets) if strong_targets else "None",
                "hub_targets_hit": ", ".join(hit_hubs) if hit_hubs else "None",
            }
        )

    out_df = pd.DataFrame(results)
    if not out_df.empty:
        out_df = out_df.sort_values(
            by=["poly_score", "strong_targets_count", "best_affinity_kcal_mol"],
            ascending=[False, False, True],
        ).reset_index(drop=True)

    return out_df
