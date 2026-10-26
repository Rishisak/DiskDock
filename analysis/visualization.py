"""Interactive Visualization Module for Step 12.
Generates responsive, interactive Vis.js network graphs and Altair enrichment charts
tailored for Streamlit with rich modern dark-theme aesthetics.
"""

from typing import List, Dict, Any, Optional
import json
import networkx as nx
import pandas as pd
import altair as alt


def generate_network_html(
    G: nx.Graph,
    title: str = "Interaction Network",
    network_type: str = "ppi",  # 'ppi' or 'target_compound'
    height_px: int = 600,
) -> str:
    """Generate a self-contained interactive Vis.js HTML string for rendering in Streamlit.

    Args:
        G: NetworkX Graph.
        title: Title of the network visualization.
        network_type: 'ppi' for Protein-Protein Interaction or 'target_compound' for Bipartite.
        height_px: Height of the visualization canvas in pixels.

    Returns:
        HTML string with embedded vis-network.min.js and interactive physics simulation.
    """
    if len(G.nodes) == 0:
        return f"""
        <div style="background:#16213e; color:#a0aec0; padding:40px; text-align:center; border-radius:12px; font-family:sans-serif;">
            <h3>No network data available</h3>
            <p>Ensure targets or compound interactions are present to render the interactome graph.</p>
        </div>
        """

    nodes_data = []
    edges_data = []

    # Calculate degrees for node sizing
    degrees = dict(G.degree())
    max_degree = max(degrees.values()) if degrees and max(degrees.values()) > 0 else 1

    for node_id in G.nodes:
        attrs = G.nodes[node_id]
        deg = degrees.get(node_id, 1)

        if network_type == "target_compound":
            ntype = attrs.get("node_type", "target")
            if ntype == "compound":
                # Emerald green for nutraceutical compounds
                color = {
                    "background": "#00e676",
                    "border": "#00b248",
                    "highlight": {"background": "#69f0ae", "border": "#ffffff"},
                }
                shape = "diamond"
                size = 18 + int((deg / max_degree) * 14)
                group = "Compound"
            else:
                # Electric cyan for target proteins
                color = {
                    "background": "#00f2fe",
                    "border": "#4facfe",
                    "highlight": {"background": "#38ef7d", "border": "#ffffff"},
                }
                shape = "dot"
                size = 14 + int((deg / max_degree) * 16)
                group = "Protein Target"
        else:
            # PPI Network: Color nodes by hub status / degree
            # High degree/hub = vibrant coral/ruby, normal = electric cyan
            is_hub = deg >= (max_degree * 0.6) and deg > 2
            if is_hub:
                color = {
                    "background": "#e94560",
                    "border": "#ff2e63",
                    "highlight": {"background": "#ff6b81", "border": "#ffffff"},
                }
                group = "Hub Protein"
                size = 22 + int((deg / max_degree) * 16)
            else:
                color = {
                    "background": "#00f2fe",
                    "border": "#00a8ff",
                    "highlight": {"background": "#70a1ff", "border": "#ffffff"},
                }
                group = "Protein Target"
                size = 14 + int((deg / max_degree) * 12)
            shape = "dot"

        nodes_data.append(
            {
                "id": str(node_id),
                "label": str(node_id),
                "title": f"<b>{node_id}</b><br>Type: {group}<br>Interactions: {deg}",
                "color": color,
                "shape": shape,
                "size": size,
                "font": {"color": "#f1f5f9", "size": 13, "face": "system-ui, -apple-system, sans-serif"},
                "borderWidth": 2,
                "shadow": {"enabled": True, "color": "rgba(0,0,0,0.4)", "size": 5},
            }
        )

    for u, v, data in G.edges(data=True):
        if network_type == "target_compound":
            aff = data.get("affinity", -7.0)
            edge_title = f"Binding Affinity: {aff:.2f} kcal/mol"
            edge_color = "rgba(0, 230, 118, 0.65)" if aff <= -8.0 else "rgba(79, 172, 254, 0.55)"
            width = max(1.5, min(5.0, (abs(aff) - 5.0) * 0.8))
        else:
            score = data.get("score", 0.4)
            edge_title = f"STRING Confidence: {score:.3f}"
            edge_color = "rgba(233, 69, 96, 0.7)" if score >= 0.7 else "rgba(112, 161, 255, 0.4)"
            width = max(1.0, min(4.5, score * 4.0))

        edges_data.append(
            {
                "from": str(u),
                "to": str(v),
                "title": edge_title,
                "color": {"color": edge_color, "highlight": "#ffffff"},
                "width": width,
                "smooth": {"type": "continuous"},
            }
        )

    nodes_json = json.dumps(nodes_data)
    edges_json = json.dumps(edges_data)

    html_code = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 0;
                background-color: #0f172a;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
                color: #e2e8f0;
                overflow: hidden;
            }}
            #network-container {{
                width: 100%;
                height: {height_px}px;
                border: 1px solid rgba(255,255,255,0.1);
                border-radius: 12px;
                background: radial-gradient(circle at 50% 50%, #1e293b 0%, #0f172a 100%);
            }}
            .hud-panel {{
                position: absolute;
                top: 14px;
                left: 14px;
                background: rgba(15, 23, 42, 0.85);
                backdrop-filter: blur(8px);
                border: 1px solid rgba(255,255,255,0.12);
                border-radius: 8px;
                padding: 8px 14px;
                font-size: 13px;
                z-index: 100;
                display: flex;
                gap: 16px;
                align-items: center;
            }}
            .legend-item {{
                display: flex;
                align-items: center;
                gap: 6px;
            }}
            .legend-dot {{
                width: 12px;
                height: 12px;
                border-radius: 50%;
                display: inline-block;
            }}
        </style>
    </head>
    <body>
        <div class="hud-panel">
            <span style="font-weight:600; color:#38bdf8;">{title}</span>
            <div class="legend-item">
                <span class="legend-dot" style="background:#00f2fe;"></span>
                <span>Targets ({len([n for n in nodes_data if n.get('shape') == 'dot'])})</span>
            </div>
            {"<div class='legend-item'><span class='legend-dot' style='background:#00e676; border-radius:2px;'></span><span>Compounds</span></div>" if network_type == "target_compound" else "<div class='legend-item'><span class='legend-dot' style='background:#e94560;'></span><span>Hub Proteins</span></div>"}
            <span style="color:#64748b;">| Drag to inspect, Scroll to zoom</span>
        </div>
        <div id="network-container"></div>

        <script type="text/javascript">
            const nodes = new vis.DataSet({nodes_json});
            const edges = new vis.DataSet({edges_json});

            const container = document.getElementById('network-container');
            const data = {{ nodes: nodes, edges: edges }};
            const options = {{
                nodes: {{
                    borderWidth: 2,
                    shadow: true
                }},
                edges: {{
                    shadow: false
                }},
                physics: {{
                    solver: 'forceAtlas2Based',
                    forceAtlas2Based: {{
                        gravitationalConstant: -55,
                        centralGravity: 0.015,
                        springLength: 100,
                        springConstant: 0.08,
                        damping: 0.6
                    }},
                    stabilization: {{
                        iterations: 120,
                        updateInterval: 25
                    }}
                }},
                interaction: {{
                    hover: true,
                    tooltipDelay: 100,
                    zoomView: true,
                    dragView: true
                }}
            }};

            const network = new vis.Network(container, data, options);
        </script>
    </body>
    </html>
    """
    return html_code


def get_enrichment_altair_chart(
    df: pd.DataFrame,
    title: str = "Enrichment Analysis",
    top_n: int = 10,
) -> Optional[alt.Chart]:
    """Generate a clean Altair horizontal bar chart showing enriched terms sorted by -log10(FDR).

    Args:
        df: Enrichment DataFrame from fetch_functional_enrichment.
        title: Title of the chart.
        top_n: Maximum number of terms to display.

    Returns:
        alt.Chart or None if DataFrame is empty.
    """
    if df.empty:
        return None

    top_df = df.head(top_n).copy()
    # Truncate overly long descriptions for clean axes
    top_df["short_desc"] = top_df["description"].apply(
        lambda x: x[:48] + "..." if len(x) > 50 else x
    )

    chart = (
        alt.Chart(top_df)
        .mark_bar(cornerRadiusTopRight=5, cornerRadiusBottomRight=5)
        .encode(
            x=alt.X(
                "neg_log10_fdr:Q",
                title="-log10(FDR)",
                scale=alt.Scale(zero=True),
            ),
            y=alt.Y(
                "short_desc:N",
                sort="-x",
                title=None,
                axis=alt.Axis(labelLimit=300),
            ),
            color=alt.Color(
                "neg_log10_fdr:Q",
                scale=alt.Scale(scheme="tealblues"),
                legend=None,
            ),
            tooltip=[
                alt.Tooltip("description:N", title="Pathway / Term"),
                alt.Tooltip("fdr:Q", format=".2e", title="FDR"),
                alt.Tooltip("neg_log10_fdr:Q", format=".2f", title="-log10(FDR)"),
                alt.Tooltip("gene_count:Q", title="Gene Count"),
                alt.Tooltip("genes:N", title="Genes"),
            ],
        )
        .properties(
            title=alt.TitleParams(text=title, fontSize=14, anchor="start"),
            height=max(220, top_n * 26),
        )
    )
    return chart
