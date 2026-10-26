import re
from pathlib import Path

import streamlit as st
import py3Dmol
import streamlit.components.v1 as components


# ============================================================
# PAGE CONFIGURATION
# ============================================================

st.set_page_config(
    page_title="DistDockNet | 3D Docking",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ============================================================
# CUSTOM CSS
# ============================================================

st.markdown(
    """
    <style>

    /* Main background */
    .stApp {
        background: #0b1020;
    }

    /* Main content */
    .block-container {
        padding-top: 2rem;
        padding-bottom: 3rem;
        max-width: 1500px;
    }

    /* Header */
    .main-title {
        font-size: 2.2rem;
        font-weight: 700;
        color: #f8fafc;
        margin-bottom: 0.2rem;
    }

    .subtitle {
        color: #94a3b8;
        font-size: 1rem;
        margin-bottom: 1.5rem;
    }

    /* Metric cards */
    .metric-card {
        background: #111827;
        border: 1px solid #263247;
        border-radius: 12px;
        padding: 18px 20px;
        min-height: 105px;
    }

    .metric-label {
        color: #94a3b8;
        font-size: 0.85rem;
        margin-bottom: 8px;
    }

    .metric-value {
        color: #f8fafc;
        font-size: 1.35rem;
        font-weight: 700;
    }

    /* Viewer container */
    .viewer-title {
        color: #f8fafc;
        font-size: 1.15rem;
        font-weight: 600;
        margin-bottom: 0.5rem;
    }

    /* Information panel */
    .info-panel {
        background: #111827;
        border: 1px solid #263247;
        border-radius: 12px;
        padding: 18px;
        margin-top: 10px;
    }

    .info-title {
        color: #f8fafc;
        font-weight: 600;
        font-size: 1rem;
        margin-bottom: 12px;
    }

    .info-row {
        display: flex;
        justify-content: space-between;
        padding: 8px 0;
        border-bottom: 1px solid #1f2937;
    }

    .info-key {
        color: #94a3b8;
    }

    .info-value {
        color: #e2e8f0;
        font-weight: 500;
    }

    /* Sidebar */
    section[data-testid="stSidebar"] {
        background: #080d1a;
        border-right: 1px solid #1e293b;
    }

    /* Buttons */
    .stButton button {
        border-radius: 8px;
    }

    /* Dataframe */
    [data-testid="stDataFrame"] {
        border-radius: 10px;
        overflow: hidden;
    }

    </style>
    """,
    unsafe_allow_html=True,
)


# ============================================================
# FILE LOCATIONS
# ============================================================

DOCKING_DIR = Path(
    r"results\ui_runs\MMP-13\Quercetin"
)

DOCKING_FILE = DOCKING_DIR / "docking_output.pdbqt"


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def find_protein_file(directory):
    """
    Find the cleaned protein PDB automatically.
    """

    preferred_names = [
        "protein_clean.pdb",
        "clean_protein.pdb",
        "clean.pdb",
        "receptor.pdb",
        "MMP-13_clean.pdb",
    ]

    for name in preferred_names:
        candidate = directory / name
        if candidate.exists():
            return candidate

    pdb_files = list(directory.glob("*.pdb"))

    if pdb_files:
        return pdb_files[0]

    return None


def read_file(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as file:
        return file.read()


def extract_vina_poses(pdbqt_text):
    """
    Extract individual MODEL blocks from Vina PDBQT output.
    """

    poses = []

    model_pattern = re.compile(
        r"MODEL\s+(\d+)(.*?)(?:ENDMDL|$)",
        re.DOTALL,
    )

    matches = model_pattern.findall(pdbqt_text)

    if matches:

        for model_number, model_text in matches:

            score = None
            rmsd_lb = None
            rmsd_ub = None

            score_match = re.search(
                r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+"
                r"(-?\d+\.\d+)\s+"
                r"(-?\d+\.\d+)",
                model_text,
            )

            if score_match:
                score = float(score_match.group(1))
                rmsd_lb = float(score_match.group(2))
                rmsd_ub = float(score_match.group(3))

            pose_text = (
                "MODEL\n"
                + model_text
                + "\nENDMDL\n"
            )

            poses.append(
                {
                    "pose": int(model_number),
                    "score": score,
                    "rmsd_lb": rmsd_lb,
                    "rmsd_ub": rmsd_ub,
                    "pdbqt": pose_text,
                }
            )

    else:

        score = None

        score_match = re.search(
            r"REMARK VINA RESULT:\s+(-?\d+\.\d+)",
            pdbqt_text,
        )

        if score_match:
            score = float(score_match.group(1))

        poses.append(
            {
                "pose": 1,
                "score": score,
                "rmsd_lb": 0.0,
                "rmsd_ub": 0.0,
                "pdbqt": pdbqt_text,
            }
        )

    return poses


def create_viewer(
    protein_data,
    ligand_data,
    protein_style="Cartoon",
    ligand_style="Sticks",
    dark_background=True,
):
    """
    Create the interactive 3Dmol viewer.
    """

    viewer = py3Dmol.view(
        width="100%",
        height=680,
    )

    # --------------------------------------------------------
    # Protein
    # --------------------------------------------------------

    viewer.addModel(
        protein_data,
        "pdb",
    )

    if protein_style == "Cartoon":

        viewer.setStyle(
            {"model": 0},
            {
                "cartoon": {
                    "color": "spectrum",
                    "opacity": 0.85,
                }
            },
        )

    elif protein_style == "Surface":

        viewer.setStyle(
            {"model": 0},
            {
                "cartoon": {
                    "color": "spectrum",
                    "opacity": 0.55,
                }
            },
        )

        viewer.addSurface(
            py3Dmol.VDW,
            {
                "opacity": 0.25,
                "color": "white",
            },
            {"model": 0},
        )

    elif protein_style == "Sticks":

        viewer.setStyle(
            {"model": 0},
            {
                "stick": {
                    "radius": 0.12,
                }
            },
        )

    # --------------------------------------------------------
    # Ligand
    # --------------------------------------------------------

    viewer.addModel(
        ligand_data,
        "pdbqt",
    )

    ligand_model = 1

    if ligand_style == "Sticks":

        viewer.setStyle(
            {"model": ligand_model},
            {
                "stick": {
                    "radius": 0.25,
                }
            },
        )

    elif ligand_style == "Ball & Stick":

        viewer.setStyle(
            {"model": ligand_model},
            {
                "stick": {
                    "radius": 0.20,
                },
                "sphere": {
                    "scale": 0.30,
                },
            },
        )

    elif ligand_style == "Spheres":

        viewer.setStyle(
            {"model": ligand_model},
            {
                "sphere": {
                    "scale": 0.45,
                }
            },
        )

    # --------------------------------------------------------
    # Background
    # --------------------------------------------------------

    if dark_background:

        viewer.setBackgroundColor(
            "#0b1020"
        )

    else:

        viewer.setBackgroundColor(
            "white"
        )

    # --------------------------------------------------------
    # Lighting / rendering
    # --------------------------------------------------------

    viewer.setProjection("perspective")

    viewer.zoomTo()

    viewer.zoom(1.15)

    viewer.spin(False)

    return viewer


# ============================================================
# HEADER
# ============================================================

st.markdown(
    '<div class="main-title">🧬 DistDockNet</div>',
    unsafe_allow_html=True,
)

st.markdown(
    '<div class="subtitle">'
    'Interactive protein–ligand docking visualization'
    '</div>',
    unsafe_allow_html=True,
)


# ============================================================
# CHECK DOCKING OUTPUT
# ============================================================

if not DOCKING_FILE.exists():

    st.error(
        "Docking output was not found."
    )

    st.code(
        str(DOCKING_FILE)
    )

    st.info(
        "Run a docking job from the DistDockNet UI first."
    )

    st.stop()


# ============================================================
# FIND PROTEIN
# ============================================================

PROTEIN_FILE = find_protein_file(
    DOCKING_DIR
)

if PROTEIN_FILE is None:

    st.error(
        "Could not find the protein PDB file."
    )

    st.write(
        "Expected a cleaned protein PDB inside:"
    )

    st.code(
        str(DOCKING_DIR)
    )

    st.stop()


# ============================================================
# LOAD DATA
# ============================================================

try:

    docking_text = read_file(
        DOCKING_FILE
    )

    protein_text = read_file(
        PROTEIN_FILE
    )

except Exception as error:

    st.error(
        f"Could not read docking files: {error}"
    )

    st.stop()


# ============================================================
# EXTRACT POSES
# ============================================================

poses = extract_vina_poses(
    docking_text
)

if not poses:

    st.error(
        "No docking poses were found."
    )

    st.stop()


# ============================================================
# SIDEBAR
# ============================================================

with st.sidebar:

    st.markdown(
        "## 🎛️ Viewer Controls"
    )

    st.markdown(
        "---"
    )

    pose_numbers = [
        pose["pose"]
        for pose in poses
    ]

    selected_pose_number = st.selectbox(
        "Docking pose",
        pose_numbers,
        index=0,
    )

    selected_pose = next(
        pose
        for pose in poses
        if pose["pose"] == selected_pose_number
    )

    st.markdown(
        "---"
    )

    protein_style = st.selectbox(
        "Protein representation",
        [
            "Cartoon",
            "Surface",
            "Sticks",
        ],
        index=0,
    )

    ligand_style = st.selectbox(
        "Ligand representation",
        [
            "Sticks",
            "Ball & Stick",
            "Spheres",
        ],
        index=0,
    )

    st.markdown(
        "---"
    )

    dark_background = st.toggle(
        "Dark viewer",
        value=True,
    )

    st.markdown(
        "---"
    )

    st.caption(
        "Drag to rotate\n\n"
        "Scroll to zoom\n\n"
        "Right-drag to move"
    )


# ============================================================
# SELECTED POSE INFORMATION
# ============================================================

score = selected_pose["score"]

if score is not None:
    score_text = f"{score:.3f} kcal/mol"
else:
    score_text = "N/A"


# ============================================================
# METRIC ROW
# ============================================================

col1, col2, col3, col4 = st.columns(4)


with col1:

    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">Protein</div>
            <div class="metric-value">MMP-13</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


with col2:

    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">Ligand</div>
            <div class="metric-value">Quercetin</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


with col3:

    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">Binding affinity</div>
            <div class="metric-value">{score_text}</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


with col4:

    st.markdown(
        f"""
        <div class="metric-card">
            <div class="metric-label">Selected pose</div>
            <div class="metric-value">
                Pose {selected_pose_number}
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )


st.write("")


# ============================================================
# MAIN LAYOUT
# ============================================================

viewer_column, information_column = st.columns(
    [3.4, 1.1],
    gap="large",
)


# ============================================================
# 3D VIEWER
# ============================================================

with viewer_column:

    st.markdown(
        '<div class="viewer-title">Interactive 3D Docking Complex</div>',
        unsafe_allow_html=True,
    )

    viewer = create_viewer(
        protein_data=protein_text,
        ligand_data=selected_pose["pdbqt"],
        protein_style=protein_style,
        ligand_style=ligand_style,
        dark_background=dark_background,
    )

    viewer_html = viewer._make_html()

    components.html(
        viewer_html,
        height=710,
        scrolling=False,
    )


# ============================================================
# INFORMATION PANEL
# ============================================================

with information_column:

    st.markdown(
        """
        <div class="info-panel">
            <div class="info-title">
                🔬 Docking information
            </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown(
        f"""
        <div class="info-row">
            <span class="info-key">Protein</span>
            <span class="info-value">MMP-13</span>
        </div>

        <div class="info-row">
            <span class="info-key">PDB</span>
            <span class="info-value">2OW9</span>
        </div>

        <div class="info-row">
            <span class="info-key">Ligand</span>
            <span class="info-value">Quercetin</span>
        </div>

        <div class="info-row">
            <span class="info-key">Pose</span>
            <span class="info-value">{selected_pose_number}</span>
        </div>

        <div class="info-row">
            <span class="info-key">Affinity</span>
            <span class="info-value">{score_text}</span>
        </div>

        <div class="info-row">
            <span class="info-key">RMSD LB</span>
            <span class="info-value">
                {selected_pose["rmsd_lb"]:.3f} Å
            </span>
        </div>

        <div class="info-row">
            <span class="info-key">RMSD UB</span>
            <span class="info-value">
                {selected_pose["rmsd_ub"]:.3f} Å
            </span>
        </div>

        <div class="info-row">
            <span class="info-key">Protein file</span>
            <span class="info-value">
                {PROTEIN_FILE.name}
            </span>
        </div>

        </div>
        """,
        unsafe_allow_html=True,
    )


# ============================================================
# POSE TABLE
# ============================================================

st.write("")

st.markdown(
    "### 📊 Docking poses"
)

pose_table = []

for pose in poses:

    pose_table.append(
        {
            "Pose": pose["pose"],
            "Binding affinity (kcal/mol)": (
                round(pose["score"], 3)
                if pose["score"] is not None
                else None
            ),
            "RMSD lower bound (Å)": (
                round(pose["rmsd_lb"], 3)
                if pose["rmsd_lb"] is not None
                else None
            ),
            "RMSD upper bound (Å)": (
                round(pose["rmsd_ub"], 3)
                if pose["rmsd_ub"] is not None
                else None
            ),
        }
    )


st.dataframe(
    pose_table,
    use_container_width=True,
    hide_index=True,
)


# ============================================================
# FOOTER
# ============================================================

st.markdown(
    """
    <div style="
        text-align:center;
        color:#64748b;
        padding-top:20px;
        font-size:0.85rem;
    ">
        DistDockNet • AutoDock Vina • 3Dmol.js
    </div>
    """,
    unsafe_allow_html=True,
)