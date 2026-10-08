# DistDockNet

### Automated Protein–Nutraceutical Molecular Docking and Analysis Pipeline

DistDockNet is an automated molecular docking pipeline designed to study interactions between disease-related protein targets and nutraceutical compounds.

The system automates the workflow from protein structure retrieval and ligand preparation to molecular docking, binding-affinity analysis, validation, and interactive 3D visualization.

---

## 📌 Overview

Molecular docking is widely used to predict how a small molecule (ligand) interacts with a target protein.

DistDockNet provides an automated workflow that integrates:

- Protein structure retrieval
- Protein preprocessing
- Compound structure generation
- Ligand preparation
- Binding-site / docking-box generation
- AutoDock Vina molecular docking
- Docking score analysis
- Pose analysis
- Redocking RMSD validation
- Interactive 3D visualization
- Target–compound ranking
- Downstream biological analysis

The goal is to reduce manual steps involved in conventional molecular docking workflows.

---

## 🧬 Project Architecture

```text
                    ┌──────────────────────┐
                    │   Protein Target     │
                    │  PDB ID / PDB File   │
                    └──────────┬───────────┘
                               │
                               ▼
                    ┌──────────────────────┐
                    │ Protein Retrieval &  │
                    │    Preprocessing     │
                    │                      │
                    │ • Structure parsing  │
                    │ • Chain selection    │
                    │ • Cleaning           │
                    │ • AltLoc handling     │
                    └──────────┬───────────┘
                               │
                               ▼
                         ┌───────────┐
                         │  Meeko    │
                         │ Receptor  │
                         │  Prep     │
                         └─────┬─────┘
                               │
                               ▼
                     Receptor PDBQT
                               │
                               │
     ┌─────────────────────────┘
     │
     │
     ▼
┌──────────────────────┐
│ Nutraceutical       │
│ Compound             │
│                      │
│ SMILES / Compound    │
│ Structure            │
└──────────┬───────────┘
           │
           ▼
┌──────────────────────┐
│ Compound Preparation │
│                      │
│ • SMILES → 3D        │
│ • H addition         │
│ • Geometry optimize  │
└──────────┬───────────┘
           │
           ▼
      ┌─────────┐
      │ Meeko   │
      │ Ligand  │
      │  Prep   │
      └────┬────┘
           │
           ▼
      Ligand PDBQT
           │
           │
           └──────────────┐
                          ▼
                 ┌─────────────────┐
                 │ Docking Box      │
                 │ Generation       │
                 │                  │
                 │ Center           │
                 │ Size             │
                 │ Search Space     │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ AutoDock Vina   │
                 │                 │
                 │ Molecular       │
                 │ Docking         │
                 └────────┬────────┘
                          │
                          ▼
                 ┌─────────────────┐
                 │ Docking Results │
                 │                 │
                 │ Affinity        │
                 │ Poses           │
                 │ RMSD            │
                 └────────┬────────┘
                          │
             ┌────────────┴────────────┐
             ▼                         ▼
   ┌──────────────────┐      ┌──────────────────┐
   │ Validation       │      │ Interaction      │
   │                  │      │ Analysis         │
   │ Redocking RMSD   │      │                  │
   │ Pose comparison  │      │ H-bonds          │
   └────────┬─────────┘      │ Hydrophobic      │
            │                │ interactions     │
            │                └────────┬─────────┘
            └────────────┬────────────┘
                         ▼
                ┌──────────────────┐
                │ Interactive 3D  │
                │ Visualization    │
                │                  │
                │ Protein + Ligand │
                │ + Docked Pose    │
                └────────┬─────────┘
                         │
                         ▼
                ┌──────────────────┐
                │ Final Ranking    │
                │                  │
                │ Protein–Compound │
                │ Interactions     │
                └──────────────────┘
