"""Comprehensive Validation Suite for Step 12: Target–Compound Intersection, PPI & Functional Analysis.
"""

import sys
from pathlib import Path
import pandas as pd

from analysis import (
    fetch_string_ppi,
    build_ppi_graph,
    calculate_network_centralities,
    normalize_target_to_gene,
    fetch_functional_enrichment,
    build_target_compound_network,
    rank_compounds_by_polypharmacology,
    generate_network_html,
    get_enrichment_altair_chart,
)

def run_step12_validation():
    print("=" * 70)
    print("         DISTDOCKNET — STEP 12 VALIDATION TEST SUITE")
    print("=" * 70)

    # ----------------------------------------------------
    # Test 1: Gene Normalization
    # ----------------------------------------------------
    raw_targets = [
        "MMP-7", "MMP-13", "MMP-19", "IL-6", "MD-2", "IL-1beta",
        "HIF-1alpha", "COL2A1", "ACAN", "ADAMTS-5", "TNF-alpha",
        "NLRP3", "IKKbeta", "SIRT1", "SOX9", "STAT3"
    ]
    normalized = [normalize_target_to_gene(t) for t in raw_targets]
    assert "TNF" in normalized, "TNF normalization failed"
    assert "IL6" in normalized, "IL6 normalization failed"
    assert "MMP13" in normalized, "MMP13 normalization failed"
    assert "IKBKB" in normalized, "IKBKB normalization failed"
    assert "LY96" in normalized, "LY96 normalization failed"
    print(f"[TEST 1/7] Gene Normalization: PASSED")
    print(f"           Cleaned 16 targets -> {normalized[:8]}...")

    # ----------------------------------------------------
    # Test 2: Live STRING DB PPI Interactome
    # ----------------------------------------------------
    edges = fetch_string_ppi(normalized, species=9606, required_score=400)
    assert len(edges) > 0, "No edges retrieved from STRING DB"
    G = build_ppi_graph(edges, user_genes=normalized)
    assert len(G.nodes) >= 16, f"Expected at least 16 nodes, got {len(G.nodes)}"
    assert len(G.edges) >= 50, f"Expected at least 50 edges, got {len(G.edges)}"
    print(f"[TEST 2/7] STRING DB PPI Interactome: PASSED")
    print(f"           Built graph with {len(G.nodes)} nodes and {len(G.edges)} verified interaction edges.")

    # ----------------------------------------------------
    # Test 3: Network Centrality & Master Hub Detection
    # ----------------------------------------------------
    centralities = calculate_network_centralities(G)
    assert not centralities.empty, "Centralities calculation returned empty DataFrame"
    assert "hub_score" in centralities.columns, "Missing hub_score"
    assert "is_hub" in centralities.columns, "Missing is_hub column"
    hubs = centralities[centralities["is_hub"]]["gene"].tolist()
    print(f"[TEST 3/7] Topological Centralities & Hub Target Detection: PASSED")
    print(f"           Master Hub Targets Detected: {', '.join(hubs)}")
    print(f"           Top 3 Hub Targets by Score:")
    for _, row in centralities.head(3).iterrows():
        print(f"             - {row['gene']}: Degree={row['degree']}, Hub Score={row['hub_score']}, Betweenness={row['betweenness_centrality']}")

    # ----------------------------------------------------
    # Test 4: Functional & Pathway Enrichment (KEGG & GO)
    # ----------------------------------------------------
    enrich = fetch_functional_enrichment(normalized, species=9606, fdr_threshold=0.05)
    kegg_df = enrich.get("KEGG Pathways", pd.DataFrame())
    gobp_df = enrich.get("GO Biological Process", pd.DataFrame())
    assert not kegg_df.empty, "KEGG enrichment empty"
    assert not gobp_df.empty, "GO-BP enrichment empty"
    assert "neg_log10_fdr" in kegg_df.columns, "Missing -log10(FDR)"
    print(f"[TEST 4/7] Functional Enrichment (KEGG & Gene Ontology): PASSED")
    print(f"           Enriched KEGG Pathways: {len(kegg_df)}")
    print(f"           Enriched GO Biological Process: {len(gobp_df)}")
    print(f"           Top KEGG Pathway: '{kegg_df.iloc[0]['description']}' (FDR = {kegg_df.iloc[0]['fdr']:.2e}, Genes: {kegg_df.iloc[0]['genes']})")

    # ----------------------------------------------------
    # Test 5: Target-Compound Intersection on Batch Data
    # ----------------------------------------------------
    batch_csv = Path("results/batch_results.csv")
    assert batch_csv.exists(), "Batch results CSV not found"
    docking_df = pd.read_csv(batch_csv)
    G_tc = build_target_compound_network(
        docking_df,
        affinity_threshold=-7.0,
        target_col="protein",
        compound_col="compound",
        affinity_col="binding_affinity"
    )
    assert len(G_tc.nodes) > 0, "Target-Compound bipartite network is empty"
    print(f"[TEST 5/7] Target-Compound Intersection Network: PASSED")
    print(f"           Constructed bipartite graph with {len(G_tc.nodes)} nodes and {len(G_tc.edges)} strong binding edges (affinity <= -7.0 kcal/mol).")

    # ----------------------------------------------------
    # Test 6: Polypharmacology Multi-Target Ranking
    # ----------------------------------------------------
    poly_df = rank_compounds_by_polypharmacology(
        docking_df,
        hub_df=centralities,
        affinity_threshold=-7.0,
        target_col="protein",
        compound_col="compound",
        affinity_col="binding_affinity"
    )
    assert not poly_df.empty, "Polypharmacology DataFrame is empty"
    print(f"[TEST 6/7] Multi-Target Polypharmacology Engine: PASSED")
    print("           Ranked Compounds:")
    for _, row in poly_df.iterrows():
        print(f"             - {row['compound']:<12}: Poly-Score={row['poly_score']:<5.1f} | Strong Targets={row['strong_targets_count']} | Best Affinity={row['best_affinity_kcal_mol']} kcal/mol | Hubs Hit={row['hub_targets_hit']}")

    # ----------------------------------------------------
    # Test 7: Interactive Visualizations (Vis.js & Altair)
    # ----------------------------------------------------
    html_ppi = generate_network_html(G, title="PPI Interactome", network_type="ppi")
    html_tc = generate_network_html(G_tc, title="Target-Compound", network_type="target_compound")
    chart = get_enrichment_altair_chart(kegg_df, title="Top KEGG")
    assert len(html_ppi) > 1000, "HTML PPI viewer generation failed"
    assert len(html_tc) > 1000, "HTML Target-Compound viewer generation failed"
    assert chart is not None, "Altair chart generation failed"
    print(f"[TEST 7/7] Interactive Visualization Generation: PASSED")
    print(f"           Generated Vis.js Canvas HTML (PPI: {len(html_ppi)} bytes, Bipartite: {len(html_tc)} bytes)")
    print(f"           Generated Altair Statistical Chart.")

    print("=" * 70)
    print("         >>> ALL 7 STEP 12 VALIDATION TESTS PASSED (100%) <<<")
    print("=" * 70)

if __name__ == "__main__":
    run_step12_validation()
