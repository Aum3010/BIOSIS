# BIOSIS: Biomolecule Interaction On Surfaces Insight System

**A Surface-Agnostic Rational Design Workbench**
*Built for the AIMS Lab Technical Challenge by Aum Pandya*

![BIOSIS Banner](https://img.shields.io/badge/Status-Prototype_Ready-success) ![Python](https://img.shields.io/badge/Python-3.10-blue) ![Streamlit](https://img.shields.io/badge/Streamlit-1.32-red)

## 🧬 Problem & Solution
Designing functionalized materials for enzyme immobilization often relies on intuition. **BIOSIS** replaces intuition with **physics-based logic**. 

It is a "Digital Twin" that bridges **Molecular Biophysics** (PDB Data) and **Surface Science** (Cheminformatics) to provide explainable insights into how proteins like *Carbonic Anhydrase* interact with any material surface.

---

## 🚀 Key Features (Meeting the Objectives)

### 1. Surface Agnostic "Surface Studio"
Unlike tools that hardcode specific materials, BIOSIS uses **SMILES strings** to parameterize *any* chemical surface.
* **Feature:** Enter `NCCN` or `C1=CC=CC=C1`.
* **Tech:** Uses **RDKit** to dynamically calculate surface descriptors (`LogP`, `TPSA`, `H-Donors`) in real-time.

### 2. The "Zinc Shield" (Active Site Protection)
I integrated the provided Zinc-distance scripts into a **Geometric Safety Engine**.
* **Feature:** A 15Å "Kill Zone" is visualized around the catalytic center.
* **Logic:** If a surface interaction is detected within this zone, the compatibility score is penalized (-40 points) to prevent enzyme deactivation.
* **Reasoning:** See "CRITICAL" warnings in the Interaction Playbook.

### 3. "The Oracle" (Literature Validation)
We don't just guess; we verify.
* **Feature:** A built-in **RAG Scout** that auto-queries OpenAlex/Crossref.
* **Logic:** Simulating "Silica"? The app runs a background search for *"Carbonic Anhydrase immobilization Silica stability"* and delivers relevant papers in the sidebar.

### 4. Explainable AI (The Playbook)
The app outputs a human-readable **Interaction Playbook**, distinguishing between:
* ✅ **Mechanistic Wins:** (e.g., "Hydrophobic Retention driven by high LogP")
* ⚠️ **Critical Risks:** (e.g., "Steric Hindrance at Active Site")

---

## 🛠️ Architecture

BIOSIS uses a modular **Functional Architecture**:

| Module | Responsibility | Tech Stack |
| :--- | :--- | :--- |
| **`bio_processor.py`** | Ingests PDBs, calculates SASA, defines Active Sites | `Biopython`, `NumPy` |
| **`surface_engine.py`** | Parameterizes surface materials | `RDKit` |
| **`interaction_lab.py`** | The "Reasoning Engine" (Scoring & Logic) | `Pandas` |
| **`literature_scout.py`** | Real-time Bibliographic Validation | `Requests`, `OpenAlex API` |

---

## 📦 How to Run

1. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
    ```
2. **Launch the Dashboard**
   ```bash
   streamlit run app.py
    ```