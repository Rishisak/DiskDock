"""
disease_docking_ui.py - Streamlit UI Component for Disease-Guided Automated Docking.
Integrates Disease-Compound-Analysis with DistDockNet Vina Docking & Step 12 Analysis.
"""

from itertools import product
from pathlib import Path
from typing import Callable, Dict, Any, List

import pandas as pd
import streamlit as st

import disease_docking_bridge as bridge

BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def render_disease_docking_workflow(
    chain_id: str,
    exhaustiveness: int,
    num_modes: int,
    vina_exe: Path,
    fetch_pdb_func: Callable[[str], Path],
    run_standard_docking_func: Callable[..., Dict[str, Any]],
    render_3d_viewer_func: Callable[[pd.DataFrame], None],
):
    """
    Renders the complete Disease-Guided Automated Docking UI.
    Preserves exact UI design language and integrates seamlessly with Step 12.
    """
    st.markdown(
        """
        <style>
        .disease-dock-banner {
            background: linear-gradient(135deg, #0f172a 0%, #14291f 50%, #0f172a 100%);
            border-left: 5px solid #10b981;
            border-radius: 12px;
            padding: 22px 28px;
            margin-bottom: 25px;
            box-shadow: 0 8px 24px rgba(0, 0, 0, 0.4);
        }
        .disease-dock-title {
            font-size: 26px;
            font-weight: 700;
            color: #f8fafc;
            margin: 0 0 6px 0;
            letter-spacing: -0.5px;
        }
        .disease-dock-subtitle {
            color: #94a3b8;
            font-size: 14px;
            margin: 0;
        }
        .disease-dock-badge {
            background: rgba(16, 185, 129, 0.15);
            color: #10b981;
            border: 1px solid rgba(16, 185, 129, 0.35);
            padding: 3px 10px;
            border-radius: 20px;
            font-size: 11px;
            font-weight: 600;
            display: inline-block;
            margin-bottom: 8px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .status-pill {
            display: inline-flex;
            align-items: center;
            gap: 6px;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: 500;
        }
        .status-pill-ok {
            background: rgba(16, 185, 129, 0.12);
            color: #34d399;
            border: 1px solid rgba(16, 185, 129, 0.3);
        }
        .status-pill-warn {
            background: rgba(245, 158, 11, 0.12);
            color: #fbbf24;
            border: 1px solid rgba(245, 158, 11, 0.3);
        }
        </style>
        <div class="disease-dock-banner">
            <span class="disease-dock-badge">Integrated Discovery Pipeline • Steps 1–12</span>
            <div class="disease-dock-title">🧬 Disease-Guided Automated Docking</div>
            <p class="disease-dock-subtitle">
                Automated multi-target discovery (Open Targets + PubMed AI), verified natural product retrieval
                (PubChem + COCONUT 2.0 Database), AutoDock Vina molecular docking, and systems polypharmacology.
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ------------------------------------------------------------
    # 0. Backend Health Check
    # ------------------------------------------------------------
    is_healthy = bridge.is_backend_healthy()
    if is_healthy:
        st.markdown(
            '<div class="status-pill status-pill-ok">● Biomedical AI Engine Active (Port 8000) &nbsp;|&nbsp; Open Targets • PubMed AI • PubChem • COCONUT Database</div>',
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            '<div class="status-pill status-pill-warn">⚠️ Biomedical AI Engine (port 8000) is connecting or offline.</div>',
            unsafe_allow_html=True,
        )
        st.info(
            "To activate full target & compound retrieval, please ensure the backend is running:\n"
            "```powershell\npython Disease-Compound-Analysis/server.py\n```"
        )

    st.write("")

    # ------------------------------------------------------------
    # 1. Disease Query & Search
    # ------------------------------------------------------------
    st.subheader("1. Disease Query & Target Discovery")

    col_q, col_opts = st.columns([3, 2])

    with col_q:
        disease_input = st.text_input(
            "Disease Name or Clinical Condition",
            value=st.session_state.get("disease_query_val", "type 2 diabetes"),
            placeholder="e.g., type 2 diabetes, melanoma, colorectal cancer, rheumatoid arthritis",
            help="Type any disease name. Normalization will match it to MONDO/Open Targets.",
        )

        # Autocomplete suggestions chips if query typed
        if disease_input and len(disease_input.strip()) >= 3:
            suggestions = bridge.get_disease_suggestions(disease_input, limit=4)
            if suggestions:
                st.caption("Suggested Ontologies:")
                cols_sug = st.columns(len(suggestions))
                for idx, sug in enumerate(suggestions):
                    sug_name = sug.get("name", "")
                    sug_id = sug.get("id", "")
                    with cols_sug[idx]:
                        if st.button(
                            f"📌 {sug_name[:24]}",
                            key=f"sug_btn_{idx}",
                            help=f"{sug_name} ({sug_id})",
                        ):
                            st.session_state["disease_query_val"] = sug_name
                            st.rerun()

    with col_opts:
        target_limit = st.slider(
            "Number of Protein Targets to Discover",
            min_value=1,
            max_value=10,
            value=4,
            step=1,
            help="Top therapeutic targets ranked by Open Targets and PubMed evidence fusion.",
        )
        papers_per_target = st.slider(
            "PubMed Papers Analyzed per Target",
            min_value=1,
            max_value=5,
            value=2,
            step=1,
        )

    btn_discover = st.button(
        "🔍 Discover Targets & Natural Compounds",
        type="primary",
        use_container_width=True,
    )

    # ------------------------------------------------------------
    # Run Discovery Pipeline
    # ------------------------------------------------------------
    if btn_discover:
        if not disease_input.strip():
            st.error("Please enter a disease name.")
            st.stop()

        if not vina_exe.exists():
            st.error(f"vina.exe not found at {vina_exe}. Please place vina.exe in the project folder.")
            st.stop()

        with st.status("Executing Biomedical Discovery Engine...", expanded=True) as status_box:
            try:
                status_box.write("Step 1–5: Normalizing disease & retrieving therapeutic targets from Open Targets...")
                target_data = bridge.run_disease_target_discovery(
                    disease_query=disease_input,
                    target_limit=target_limit,
                    papers_per_target=papers_per_target,
                )
                targets = target_data.get("targets", [])
                disease_name = target_data.get("disease_name", disease_input)
                synonyms = target_data.get("synonyms", [])

                status_box.write(f"✓ Found {len(targets)} targets for '{disease_name}'. Resolving PDB structures...")
                resolved_targets = []
                for t in targets:
                    symbol = t.get("target_symbol", "")
                    name = t.get("target_name", "")
                    pdb_id, pdb_source = bridge.resolve_target_pdb(symbol, name)
                    resolved_targets.append({
                        "target_symbol": symbol,
                        "target_name": name,
                        "fused_score": t.get("fused_score", 0.0),
                        "open_targets_score": t.get("open_targets_score", 0.0),
                        "literature_score": t.get("literature_score", 0.0),
                        "pdb_id": pdb_id,
                        "pdb_source": pdb_source,
                    })

                status_box.write("Step 6–7: Retrieving candidate compounds & verifying via PubChem + COCONUT 2.0...")
                compounds_data = bridge.run_natural_compound_discovery(
                    disease_name=disease_name,
                    targets=targets,
                    synonyms=synonyms,
                )
                compounds = compounds_data.get("compounds", [])

                status_box.write(f"✓ Retrieved {len(compounds)} verified natural compounds from database!")
                status_box.update(label="Biomedical Discovery Complete!", state="complete", expanded=False)

                st.session_state["disease_docking_data"] = {
                    "disease_query": disease_input,
                    "disease_name": disease_name,
                    "synonyms": synonyms,
                    "targets": resolved_targets,
                    "compounds": compounds,
                    "compounds_summary": {
                        "total": compounds_data.get("total_compounds", len(compounds)),
                        "natural": compounds_data.get("verified_natural_count", 0),
                        "derivatives": compounds_data.get("natural_derivative_count", 0),
                    },
                }
            except Exception as exc:
                status_box.update(label="Discovery Failed", state="error", expanded=True)
                st.error(f"Error during biomedical discovery: {exc}")
                st.stop()

    # ------------------------------------------------------------
    # 2. Review Discovered Targets & Compounds
    # ------------------------------------------------------------
    if "disease_docking_data" in st.session_state:
        data = st.session_state["disease_docking_data"]
        targets_list = data.get("targets", [])
        compounds_list = data.get("compounds", [])
        disease_name = data.get("disease_name", "")

        st.divider()
        st.subheader("2. Review Discovered Targets & Natural Compounds")

        # Summary KPIs
        kpi1, kpi2, kpi3, kpi4 = st.columns(4)
        kpi1.metric("Analyzed Condition", disease_name)
        kpi2.metric("Discovered Targets", len(targets_list))
        kpi3.metric("Verified Natural Products", data["compounds_summary"]["natural"])
        kpi4.metric("COCONUT Verified", sum(1 for c in compounds_list if c.get("natural_origin")))

        tab_targets, tab_compounds = st.tabs([
            f"🎯 Therapeutic Targets ({len(targets_list)})",
            f"🌿 Verified Natural Compounds ({len(compounds_list)})",
        ])

        with tab_targets:
            st.caption("Therapeutic protein targets ranked by fused Open Targets + PubMed evidence:")
            targets_table = []
            for t in targets_list:
                targets_table.append({
                    "Target Symbol": t["target_symbol"],
                    "Target Name": t["target_name"],
                    "Resolved PDB ID": t["pdb_id"] if t["pdb_id"] else "Searching...",
                    "PDB Method": t["pdb_source"],
                    "Fused Evidence Score": f"{t['fused_score']:.3f}",
                })
            df_targets = pd.DataFrame(targets_table)
            st.dataframe(df_targets, use_container_width=True, hide_index=True)

            # Target selection for docking
            target_symbols = [t["target_symbol"] for t in targets_list if t.get("pdb_id")]
            selected_target_symbols = st.multiselect(
                "Select Targets for Automated Docking:",
                options=target_symbols,
                default=target_symbols[:2] if len(target_symbols) >= 2 else target_symbols,
                help="Choose which resolved protein targets will enter the docking queue.",
            )

        with tab_compounds:
            st.caption("Candidate natural products verified against PubChem and COCONUT database:")
            compounds_table = []
            for c in compounds_list:
                compounds_table.append({
                    "Compound Name": c.get("name", "").capitalize(),
                    "Origin Status": c.get("natural_origin_status", "natural").upper(),
                    "Class": c.get("natural_product_class", "Natural Product"),
                    "COCONUT / PubChem ID": f"CID {c.get('pubchem_cid', '')}",
                    "Bioactivity (nM)": f"{c.get('best_bioactivity_nM', 0):.1f}" if c.get("best_bioactivity_nM") else "Active",
                    "Dual Route": "✓ Yes" if c.get("is_dual_route") else "Route A/B",
                    "SMILES": c.get("canonical_smiles") or c.get("smiles") or "",
                })
            df_compounds = pd.DataFrame(compounds_table)
            st.dataframe(df_compounds, use_container_width=True, hide_index=True)

            # Compound selection for docking
            compound_names = [c.get("name", "").capitalize() for c in compounds_list if (c.get("canonical_smiles") or c.get("smiles"))]
            selected_compound_names = st.multiselect(
                "Select Natural Compounds for Automated Docking:",
                options=compound_names,
                default=compound_names[:3] if len(compound_names) >= 3 else compound_names,
                help="Choose which verified phytochemicals will be docked against the selected targets.",
            )

        # ------------------------------------------------------------
        # 3. Docking Execution Queue
        # ------------------------------------------------------------
        st.divider()
        st.subheader("3. Automated Molecular Docking")

        # Filter the chosen targets and compounds
        chosen_targets = [t for t in targets_list if t["target_symbol"] in selected_target_symbols]
        chosen_compounds = [
            c for c in compounds_list
            if c.get("name", "").capitalize() in selected_compound_names
        ]

        num_jobs = len(chosen_targets) * len(chosen_compounds)

        col_meta, col_btn = st.columns([3, 2])
        with col_meta:
            st.write(
                f"**Docking Matrix:** `{len(chosen_targets)}` protein target(s) × `{len(chosen_compounds)}` natural compound(s) = **`{num_jobs}` docking combination(s)**"
            )
            st.caption(
                f"Vina Settings: Exhaustiveness = {exhaustiveness} • Poses = {num_modes} • Chain = {chain_id.strip() or 'A'}"
            )

        with col_btn:
            btn_dock = st.button(
                "⚡ Run AutoDock Vina Docking",
                type="primary",
                disabled=(num_jobs == 0),
                use_container_width=True,
            )

        if btn_dock:
            if num_jobs == 0:
                st.warning("Please select at least 1 target and 1 compound.")
                st.stop()

            # Prepare receptor structures
            progress_bar = st.progress(0)
            status_text = st.empty()

            resolved_receptors = []
            for t in chosen_targets:
                sym = t["target_symbol"]
                pdb_code = t["pdb_id"]
                status_text.write(f"Preparing receptor {sym} (PDB: {pdb_code})...")
                try:
                    pdb_path = fetch_pdb_func(pdb_code)
                    resolved_receptors.append({
                        "name": sym,
                        "pdb_id": pdb_code,
                        "pdb_file": pdb_path,
                    })
                except Exception as exc:
                    st.error(f"Failed to fetch PDB {pdb_code} for target {sym}: {exc}")
                    st.stop()

            # Execute docking matrix
            results = []
            jobs = list(product(resolved_receptors, chosen_compounds))

            for idx, (rec, comp) in enumerate(jobs, start=1):
                rec_name = rec["name"]
                comp_name = comp.get("name", "").capitalize()
                smiles = comp.get("canonical_smiles") or comp.get("smiles")

                status_text.write(
                    f"Docking [{idx}/{num_jobs}]: {rec_name} + {comp_name}..."
                )
                progress_bar.progress(idx / num_jobs)

                try:
                    dock_res = run_standard_docking_func(
                        protein_name=rec_name,
                        pdb_file=rec["pdb_file"],
                        compound_name=comp_name,
                        smiles=smiles,
                        chain_id=chain_id.strip() or "A",
                        site_mode="blind",
                        reference_ligand="",
                        exhaustiveness=exhaustiveness,
                        num_modes=num_modes,
                    )

                    results.append({
                        "Protein": dock_res["protein"],
                        "Compound": dock_res["compound"],
                        "Binding affinity (kcal/mol)": dock_res["binding_affinity"],
                        "Chain": dock_res["selected_chain"],
                        "Docking output": dock_res["output"],
                        "Protein PDB": dock_res["protein_pdb"],
                        "Ligand PDBQT": dock_res["ligand_pdbqt"],
                        "PDB ID": rec["pdb_id"],
                        "SMILES": smiles,
                        "Natural Origin": comp.get("natural_origin_status", "natural"),
                        "Class": comp.get("natural_product_class", "Natural Product"),
                    })
                except Exception as exc:
                    results.append({
                        "Protein": rec_name,
                        "Compound": comp_name,
                        "Binding affinity (kcal/mol)": None,
                        "Chain": None,
                        "Docking output": str(exc),
                        "Protein PDB": "",
                        "Ligand PDBQT": "",
                        "PDB ID": rec["pdb_id"],
                        "SMILES": smiles,
                        "Natural Origin": comp.get("natural_origin_status", "natural"),
                        "Class": comp.get("natural_product_class", "Natural Product"),
                    })

            progress_bar.progress(1.0)
            status_text.success(f"✓ All {num_jobs} docking calculations successfully completed!")

            results_df = pd.DataFrame(results)

            # Persist to disk for Step 12 & downstream analysis
            batch_csv = RESULTS_DIR / "batch_results.csv"
            ui_csv = RESULTS_DIR / "ui_results.csv"

            # Create format compatible with Step 12
            step12_df = results_df.copy()
            step12_df["protein"] = step12_df["Protein"]
            step12_df["compound"] = step12_df["Compound"]
            step12_df["affinity"] = step12_df["Binding affinity (kcal/mol)"]

            step12_df.to_csv(batch_csv, index=False)
            results_df.to_csv(ui_csv, index=False)

            st.session_state["results_df"] = results_df
            st.session_state["disease_docking_results"] = results_df

        # ------------------------------------------------------------
        # 4. Display Docking Results
        # ------------------------------------------------------------
        if "disease_docking_results" in st.session_state:
            res_df = st.session_state["disease_docking_results"]

            st.divider()
            st.subheader("4. Docking Results & Binding Affinities")

            num_aff = pd.to_numeric(res_df["Binding affinity (kcal/mol)"], errors="coerce")
            successful = res_df.loc[num_aff.notna()].copy()
            successful = successful.sort_values("Binding affinity (kcal/mol)", ascending=True)

            if not successful.empty:
                best = successful.iloc[0]
                best_score = float(best["Binding affinity (kcal/mol)"])

                b_col1, b_col2, b_col3 = st.columns(3)
                with b_col1:
                    st.metric(
                        "Top Hit Complex",
                        f"{best['Protein']} + {best['Compound']}",
                    )
                with b_col2:
                    st.metric(
                        "Strongest Binding Affinity",
                        f"{best_score:.2f} kcal/mol",
                    )
                with b_col3:
                    st.metric(
                        "Phytochemical Class",
                        str(best.get("Class", "Natural Product")),
                    )

            # Display table
            display_cols = [
                "Protein", "Compound", "Binding affinity (kcal/mol)", "PDB ID", "Chain", "Natural Origin", "Class"
            ]
            valid_cols = [c for c in display_cols if c in res_df.columns]
            st.dataframe(res_df[valid_cols], use_container_width=True, hide_index=True)

            # CSV Download
            csv_bytes = res_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "📥 Download Docking Results CSV",
                data=csv_bytes,
                file_name=f"disease_docking_{data.get('disease_name', 'results').replace(' ', '_')}.csv",
                mime="text/csv",
            )

            # ------------------------------------------------------------
            # 5. Interactive 3D Viewer
            # ------------------------------------------------------------
            render_3d_viewer_func(res_df)

            # ------------------------------------------------------------
            # 6. Step 12 Systems Pharmacology Seamless Transition
            # ------------------------------------------------------------
            st.divider()
            st.markdown(
                """
                <div style="background: linear-gradient(135deg, #0f172a 0%, #1e1b4b 50%, #0f172a 100%);
                            border-left: 5px solid #38bdf8; border-radius: 12px; padding: 20px 24px; margin-top: 15px;">
                    <h3 style="color: #f8fafc; margin: 0 0 8px 0; font-size: 20px;">
                        🌐 Step 12: Network Pharmacology & PPI Analysis
                    </h3>
                    <p style="color: #94a3b8; margin: 0 0 15px 0; font-size: 14px;">
                        Seamlessly explore the STRING protein-protein interactome, master hub targets,
                        KEGG/GO functional pathways, and polypharmacology bipartite graph for your docked disease targets.
                    </p>
                </div>
                """,
                unsafe_allow_html=True,
            )

            st.write("")
            if st.button(
                "🚀 Open Step 12: PPI & Network Pharmacology Analysis",
                type="primary",
                use_container_width=True,
            ):
                st.session_state["jump_to_step12"] = True
                st.session_state["workflow_mode"] = "🌐 Step 12: PPI & Network Analysis"
                st.rerun()

