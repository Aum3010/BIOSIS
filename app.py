import os
import time
import hashlib
import json
import datetime
import altair as alt
import pandas as pd
import streamlit as st
from dotenv import load_dotenv
from modules.bio_processor import process_protein_structure, process_pdb_string, mutate_residue, get_rotated_pdb, map_metric_to_bfactor
from modules.surface_engine import SurfaceEngine
from modules.interaction_lab import InteractionLab
from modules.literature_scout import LiteratureScout
from modules.visualizer import SurfaceVisualizer
from modules.optimizer import GeneticOptimizer
from modules.generator import MaterialGenerator
from stmol import showmol
import py3Dmol

load_dotenv()

def draw_help_center():
    """Renders the pragmatic Help Guide mapped to UI locations."""
    with st.popover("Help & Guide", use_container_width=True):
        st.markdown("### 📘 BIOSIS User Manual")
        st.caption("Comprehensive guide to buttons, physics, and workflows.")
        
        tab_side, tab_view, tab_data, tab_tools, tab_file = st.tabs([
            "🎛️ Sidebar", 
            "🧊 3D Viewer", 
            "📈 Graphs & Metrics", 
            "🧠 AI & Tuning", 
            "📂 Downloads"
        ])
        
        with tab_side:
            st.markdown("#### 1. Sidebar Configuration")
            st.info("Controls for the Physical Reactor Conditions.")
            
            st.markdown("##### **A. Material Selection**")
            st.markdown("""
            * **What it is:** Defines the chemical properties (Charge, Hydrophobicity) of the sensor surface.
            * **Purpose:** Matching the protein to the surface prevents denaturation (unfolding).
            * **Standard Options:**
                * `Gold`
                * `Silica`
                * `Graphene`
            * **Advanced Options:**
                * `Custom`: Enter any **SMILES** string (e.g., `NCC` for Ethylamine) to simulate your specific proprietary material.
                * `AI Generated`: An **Inverse Design Engine**. Select a goal (e.g., *Hydrophobic*, *Hydrophilic*, or *Positive*), and the AI will find a new chemical scaffold optimized for that property using genetic algorithm.
            """)
            
            st.divider()
            
            st.markdown("##### **B. Linker Strategy (Bioconjugation)**")
            st.markdown("""
            * **What it is:** A chemical tether that lifts the protein off the surface floor.
            * **Purpose:**  Solves the "Landing Problem." Proteins touching the floor often lose activity due to steric hindrance or denaturation.
            * **The Trade-Off:**
                * **Short (0-15Å):** High control over orientation, but higher risk of surface impact.
                * **Long (40Å+):** High safety (protein floats freely), but high "Entropy" (wobbles like a balloon), reducing sensor precision.
            * **Action:** Start with `PEG-4`. If score is low due to "Steric", switch to `PEG-12`.
            """)

            st.divider()

            st.markdown("##### **C. Crowding & Reactor**")
            st.markdown("""
            * **Target Density:** The concentration of protein on the chip (pmol/cm²).
            * **Purpose:**  Real sensors are packed monolayers. If packed too tight, proteins block each other's active sites ("Steric Choking").
            * **pH Slider:** Changes the protonation state of amino acids (His, Lys, Asp, Glu).
            * **Example:** *At pH 5, Histidine is Positive (+). At pH 7, it is Neutral (0). This drastically changes binding.*
            """)
            
        with tab_view:
            st.markdown("#### 2. The 3D Viewer (Left Panel)")
            st.info("Visualizing the 'Digital Twin' simulation.")
            
            st.markdown("##### **A. Single (Chemistry)**")
            st.markdown("""
            * **What it shows:** Standard atomic representation (CPK Coloring).
            * **Purpose:** Verify physical orientation.
            * **Key Detail:** Look for **Red Spheres**. These represent Zinc atoms (Active Sites). They must point *up* (away from surface) for the sensor to work.
            """)

            st.markdown("##### **B. Single (Evolution)**")
            st.markdown("""
            * **What it shows:** Evolutionary Conservation scores mapped to the structure.
            * **Purpose:** Identifies which parts of the protein are "vital" vs. "modifiable."
            * **Color Code:**
                * 🔴 **Red (Conserved):**  Functional core. **DO NOT MUTATE.**
                * 🔵 **Blue (Variable):** Structural scaffold. Safe to engineer/mutate.
            """)

            st.markdown("##### **C. Single (Physics)**")
            st.markdown("""
            * **What it shows:** Eisenberg Hydrophobicity Scale.
            * **Purpose:** detecting "Hydrophobic Mismatch."
            * **Color Code:**
                * 🔴 **Red:** Polar (Water-loving).
                * 🔵 **Blue:** Hydrophobic (Oily/Sticky).
            """)

            st.markdown("##### **D. Lattice (Crowding)**")
            st.markdown("""
            * **What it shows:** A 3x3 grid of protein clones.
            * **Purpose:** Visualizing lateral steric hindrance.
            * **Action:** Drag the `Density` slider. If the clones intersect, you are over-packing the sensor.
            """)

        with tab_data:
            st.markdown("#### 3. Graphs & Scoreboard (Right Panel)")
            st.info("Interpreting the simulation metrics.")
            
            st.markdown("##### **A. The Scoreboard**")
            st.markdown("""
            * **BIOSIS Score (0-100):** A weighted probability of sensor success.
                * *<50:* Likely failure (Denaturation or bad orientation).
                * *>80:* Manufacturable candidate.
            * **Surface Hydro:** LogP value. High = Oily. Low = Watery.
            * **Cons. Risk:** "High" means the interface involves evolutionarily conserved residues.
            """)

            st.markdown("##### **B. Literature Scout**")
            st.markdown("""
            * **What it is:** A real-time search engine connected to OpenAlex/Crossref.
            * **Purpose:** "Meant to find specific prior research on similar protein-surface systems."
            * **How it works:** It constructs a query like `"[Protein] + [Surface] + Immobilization"` to find relevant prior research.
            """)

            st.markdown("##### **C. Stability Profile**")
            st.markdown("""
            * **What it is:** A digital pH titration curve.
            * **Purpose:** Determining the robust operating range.
            * **Interpretation:**
                * **Flat Line:** The sensor works across all pH levels.
                * **Sharp Drop:** The sensor will fail if pH shifts slightly.
            """)

        with tab_tools:
            st.markdown("#### 4. AI & Tuning (Right Tabs)")
            st.info("Advanced tools for fixing and refining the model.")
            
            st.markdown("##### **A. AI Architect (Mutagenesis)**")
            st.markdown("""
            * **What it is:** A Generative Design agent powered by GPT-4o.
            * **Purpose:**  Solving physics problems that parameters can't fix.
            * **Trigger:** If the physics engine detects a clash (e.g., "Electrostatic Repulsion").
            * **Action:** The AI suggests a specific mutation (e.g., `ASP-123 -> LYS`) to flip the charge.
            * **Safety Gate:** The AI is strictly blocked from suggesting mutations on Conserved (Red) residues.
            """)

            st.markdown("##### **B. Lab-in-the-Loop (Calibration)**")
            st.markdown("""
            * **What it is:** A Bayesian Weight Refinement engine.
            * **Purpose:** Customizing the physics engine to *your* specific lab reality.
            * **Scenario:** The app predicts 80% score, but your lab result is 20%.
            * **Workflow:**
                1.  Enter "20" in the **Actual Lab Yield** box.
                2.  Click **Add to Training Set**.
                3.  Click **Run Refinement**.
            * **Result:** The math engine "re-weights" the importance of Sterics vs. Electrostatics to minimize the error for your future runs.
            """)

        with tab_file:
            st.markdown("#### 5. Downloads (Export Tab)")
            st.info("Translating code into physical reality.")
            
            st.markdown("##### **A. Compliance Report**")
            st.markdown("""
            * **What it is:** A read-only audit trail.
            * **Purpose:** Regulatory documentation (ISO/FDA).
            * **Security:** Includes a **SHA-256 Hash**. Any alteration to the file breaks the signature, proving data integrity.
            """)

            st.markdown("##### **B. Robot Protocol**")
            st.markdown("""
            * **What it is:** A Python script (`.py`).
            * **Purpose:**  Automating the wet-lab work.
            * **Compatibility:** Native support for **Opentrons OT-2** robots.
            * **Logic:** It reads your `Density` settings to calculate exactly how many microliters of protein/linker to pipette.
            """)
            
def draw_landing_page():
    """Renders the detailed, product-grade onboarding dashboard."""        
    col1, col2, col3 = st.columns(3)
    with col1:
        with st.container(border=True):
            st.markdown("### 🧬 Phase 1: Design")
            st.markdown("""
            **Simulate the Thermodynamics:**
            * **Physics Engine:** Calculates Electrostatic, Hydrophobic, and Steric interaction energies.
            * **Dynamic Reactor:**  Adjust pH (4-10) to see how protein charge (pI) affects binding.
            * **Safety Checks:** Automatically detects if the Active Site (Zinc) is blocked by the surface.
            
            👉 *Try: Input PDB `1V9E` (Streptavidin) and test it on `Gold` vs `Silica` surfaces.*
            """)
    
    with col2:
        with st.container(border=True):
            st.markdown("### 🏗️ Phase 2: Engineer")
            st.markdown("""
            * **Linker Lab:** Don't let proteins hit the floor. Add `PEG-12` tethers to prevent denaturation.
            * **Crowding Simulator:**  Calculate the maximum `pmol/cm²` density before proteins block each other ("Steric Choking").
            
            👉 *Try: Go to `Crowding` tab and set Density to `10 pmol/cm²` to see a lattice clash.*
            """)
        
    with col3:
        with st.container(border=True):
            st.markdown("### 🏭 Phase 3: Build")
            st.markdown("""
            **Manufacturing & Tuning:**
            * **Literature Scout:** Find specific prior research on similar protein-surface systems.
            * **AI Architect:** GPT-4o suggests point mutations (e.g., `ASP->LYS`) to fix binding faults.
            * **Lab-in-the-Loop:** Feed real lab results back into the model to calibrate its accuracy.
            * **Robot Export:** Download valid **Opentrons** scripts to automate the wet lab.
            
            👉 *Try: Go to `Export` to get your SHA-256 signed Compliance Report.*
            """)

    st.markdown("---")
    st.success("👈 **Get Started:** Configure your Protein (PDB) and Surface Material in the Sidebar to begin the simulation.")

def generate_audit_report(pdb, surface, linker, ph, score, warnings):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    status = "APPROVED" if score > 70 else "FLAGGED"
    content = f"""
===================================================
BIOSIS ARCHITECT - COMPLIANCE AUDIT REPORT
===================================================
Date: {timestamp}
Session ID: {hashlib.md5(timestamp.encode()).hexdigest()[:8]}
Operator: Session_User

1. SYSTEM CONFIGURATION
-----------------------
Protein Target:  {pdb}
Surface Material: {surface}
Bioconjugation:   {linker}
Reactor State:    pH {ph}

2. SIMULATION METRICS
---------------------
BIOSIS Score:     {score:.1f}/100
Risk Assessment:  {status}
Critical Flags:   {len(warnings)}
{chr(10).join([f" - [WARN] {w}" for w in warnings])}

3. INTEGRITY SIGNATURE
----------------------
SHA-256: {hashlib.sha256(f"{pdb}{surface}{score}".encode()).hexdigest()}
===================================================
"""
    return content

def generate_opentrons_protocol(pdb, linker_type, ph, vol_ul=100):
    return f"""
from opentrons import protocol_api

metadata = {{
    'protocolName': 'BIOSIS Immobilization: {pdb}',
    'description': 'Protocol for {linker_type} at pH {ph}',
    'apiLevel': '2.13'
}}

def run(protocol: protocol_api.ProtocolContext):
    plate = protocol.load_labware('corning_96_wellplate_360ul_flat', 1)
    tiprack = protocol.load_labware('opentrons_96_tiprack_300ul', 2)
    p300 = protocol.load_instrument('p300_single_gen2', 'right', tip_racks=[tiprack])

    protein_stock = plate['A1']
    buffer_stock = plate['A2']
    linker_stock = plate['A3']
    targets = plate.rows()[1][:6] 

    protocol.comment('Adding Buffer...')
    p300.distribute(50, buffer_stock, targets)
    
    if '{linker_type}' != 'Direct Adsorption':
        protocol.comment('Adding Linker...')
        p300.transfer(10, linker_stock, targets, mix_after=(3, 20))
        protocol.delay(minutes=15)

    protocol.comment('Adding Protein...')
    p300.transfer({vol_ul/2}, protein_stock, targets, mix_after=(5, 50))
"""

st.set_page_config(page_title="BIOSIS", layout="wide", page_icon="🧬")

@st.cache_resource
def get_static_engines():
    return SurfaceEngine(), InteractionLab(), LiteratureScout(), SurfaceVisualizer(), MaterialGenerator()

chem_engine, lab_brain, scout, vis_engine, generator = get_static_engines()

if 'current_pdb' not in st.session_state: st.session_state['current_pdb'] = "1V9E"
if 'mutated_pdb_str' not in st.session_state: st.session_state['mutated_pdb_str'] = None
if 'generated_surface' not in st.session_state: st.session_state['generated_surface'] = None
if 'calibration_data' not in st.session_state: st.session_state['calibration_data'] = []

# --- SIDEBAR ---
with st.sidebar:
    api_key = os.getenv("OPENAI_API_KEY")    
    st.subheader("Protein")
    pdb_input = st.text_input("PDB ID", st.session_state['current_pdb']).upper()
    if pdb_input != st.session_state['current_pdb']:
        st.session_state['current_pdb'] = pdb_input
        st.session_state['mutated_pdb_str'] = None
        st.rerun()

    if st.session_state['mutated_pdb_str']:
        st.warning("Variant Active")
        if st.button("Reset to Wild Type"):
            st.session_state['mutated_pdb_str'] = None
            st.rerun()

    st.subheader("Surface Chemistry")
    s_choice = st.selectbox("Material", ["Graphene", "Silica", "Gold", "Custom", "AI Generated"])
    
    if s_choice == "AI Generated":
        with st.expander("Discovery Engine", expanded=True):
            goal = st.selectbox("Optimization Goal", ["hydrophobic", "hydrophilic", "positive"])
            scaffold = st.text_input("Base Scaffold", "c1ccccc1")
            if st.button("Screen Library"):
                with st.spinner("Screening..."):
                    cands = generator.generate_candidates(scaffold, goal)
                    st.session_state['generated_surface'] = cands[0]['smiles']
                    st.success(f"Selected: {cands[0]['name']}")
                    time.sleep(1)
                    st.rerun()
    
    base_smiles = "c1ccccc1"
    if s_choice == "Silica": base_smiles = "O=[Si]=O"
    elif s_choice == "Gold": base_smiles = "[Au]"
    elif s_choice == "Custom": base_smiles = st.text_input("SMILES", "NCC")
    elif s_choice == "AI Generated" and st.session_state['generated_surface']:
        base_smiles = st.session_state['generated_surface']
        st.caption(f"SMILES: {base_smiles[:20]}...")

    st.markdown("**Linker Lab**")
    linker_label = st.selectbox("Bioconjugation Strategy", 
                                ["Direct Adsorption (0 Å)", 
                                 "PEG-4 Spacer (Flexible, 15 Å)", 
                                 "PEG-12 Tether (Flexible, 40 Å)", 
                                 "Alkane Chain (Rigid, 20 Å)"])
    
    linker_len = 0.0
    linker_flex = "flexible"
    if "PEG-4" in linker_label: linker_len = 15.0
    elif "PEG-12" in linker_label: linker_len = 40.0
    elif "Alkane" in linker_label: 
        linker_len = 20.0
        linker_flex = "rigid"

    st.subheader("Reactor")
    ph = st.slider("pH", 4.0, 10.0, 7.4)
    density = st.select_slider("Crowding", ["Single", "High (Monolayer)"])
    topology = st.selectbox("Topology", ["Flat", "Nanopore", "Nanopillar"])
    
    with st.expander("Crowding Settings"):
        target_density = st.slider("Density (pmol/cm²)", 0.1, 20.0, 5.0)
    
    st.markdown("---")    
    draw_help_center()
    
    run_sim = st.button("🚀 Run Simulation", type="primary", use_container_width=True)

st.title("BIOSIS (Biomolecule Interaction On Surfaces Insight System)")

if not run_sim and not st.session_state.get('last_sim'):
    draw_landing_page()
else:
    st.session_state['last_sim'] = True
    
    with st.spinner("Processing ..."):
        if st.session_state['mutated_pdb_str']:
            df_surface, pdb_str, zinc_coords = process_pdb_string(st.session_state['mutated_pdb_str'])
            meta = {"title": "Engineered Variant", "method": "In-Silico"}
        else:
            df_surface, pdb_str, zinc_coords, meta = process_protein_structure(pdb_input)

        surface_props = chem_engine.process_smiles(base_smiles)
        optimizer = GeneticOptimizer(df_surface, zinc_coords)
        t_type = "positive" if surface_props['tpsa'] > 40 else "hydrophobic"
        best_rot, _, _ = optimizer.run(surface_type=t_type)
        pdb_rotated = get_rotated_pdb(pdb_str, best_rot)  
        score, reasons, warnings, ai_suggestions, components = lab_brain.calculate_interface_score(
            df_surface, surface_props, best_rot, ph, density, topology, api_key,
            linker_length=linker_len, linker_type=linker_flex
        )

    st.caption("Legend: 🤖 Inferred using AI | 🧮 Calculated | 📏 Measured")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("BIOSIS Score 🤖", f"{score:.0f}", help="Success Probability")
    c2.metric("Surface Hydro. 🧮", f"{surface_props.get('logP',0):.1f}")
    c3.metric("Cons. Risk 🧮", "High" if df_surface['Conservation'].mean() > 0.8 else "Low")
    c4.metric("Linker 📏", f"{linker_len}Å", delta=linker_flex)

    st.markdown("---")

    col_left, col_right = st.columns([1.8, 1.2])

    with col_left:
        c_v1, c_v2 = st.columns([3, 2])
        with c_v1: st.subheader("3D Reality Simulator")
        with c_v2: 
            view_mode = st.radio("View Mode:", 
                                ["Single (Chemistry)", "Single (Evolution)", "Single (Physics)", "Lattice (Crowding)"], 
                                horizontal=True, label_visibility="collapsed")

        view = py3Dmol.view(height=600)

        if "Lattice" in view_mode:
            spacing_val, _, _, _ = lab_brain.calculate_crowding_metrics(df_surface, target_density)
            lattice_pdb = vis_engine.generate_protein_lattice(pdb_rotated, spacing_val)
            if lattice_pdb:
                view.addModel(lattice_pdb, 'pdb')
                view.setStyle({'chain': 'A'}, {'cartoon': {'color': 'green'}})
                view.setStyle({'chain': 'B'}, {'cartoon': {'color': 'white', 'opacity': 0.5}})
                st.caption(f"Viewing 3x3 Lattice at {target_density} pmol/cm² (Spacing: {spacing_val:.1f}Å)")
        else:
            viz_pdb = pdb_rotated
            if "Evolution" in view_mode:
                viz_pdb = map_metric_to_bfactor(pdb_rotated, df_surface, "Conservation")
            elif "Physics" in view_mode:
                viz_pdb = map_metric_to_bfactor(pdb_rotated, df_surface, "Hydrophobicity")

            view.addModel(viz_pdb, 'pdb')
            
            if "Chemistry" in view_mode:
                view.setStyle({'cartoon': {'color': 'white'}})
                for z in zinc_coords:
                     view.addSphere({'center':{'x':z[0],'y':z[1],'z':z[2]}, 'radius': 2.0, 'color': 'red'})
                view.addStyle({'resn': 'LYS'}, {'stick': {'colorscheme': 'blueCarbon'}})
            elif "Evolution" in view_mode:
                view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'rwb', 'min': 0, 'max': 100}}})
                st.caption("🔴 Red = Conserved | 🔵 Blue = Variable")
            elif "Physics" in view_mode:
                view.setStyle({'cartoon': {'colorscheme': {'prop': 'b', 'gradient': 'rwb', 'min': 0, 'max': 100}}})
                st.caption("🔴 Red = Polar | 🔵 Blue = Hydrophobic")

        surface_offset = -15.0 - linker_len
        surface_pdb = vis_engine.generate_topology(base_smiles, shape=topology, density=density, z_offset=surface_offset)
        if surface_pdb:
            view.addModel(surface_pdb, 'pdb')
            view.setStyle({'model': -1}, {'stick': {}})
            
        if linker_len > 0:
            tether_color = 'red' if linker_flex == "flexible" else 'yellow'
            style = {'dashed': True} if linker_flex == "flexible" else {}
            view.addCylinder({
                'start': {'x': 0, 'y': 0, 'z': -8}, 
                'end': {'x': 0, 'y': 0, 'z': surface_offset},
                'radius': 0.5, 'color': tether_color, **style
            })
            
        view.zoomTo()
        showmol(view, height=600)

    with col_right:
        tab_analysis, tab_crowd, tab_actions, tab_calib, tab_export = st.tabs(["🧠 Analysis", "👥 Crowding", "🛠️ Actions", "⚙️ Lab Loop", "📂 Export"])
        
        with tab_analysis:
            scout.render_feed(meta.get('title', pdb_input), s_choice)
            st.markdown("#### Stability Profile")
            res = [{"pH": p, "Score": lab_brain.calculate_interface_score(df_surface, surface_props, best_rot, float(p), density, topology, api_key, linker_length=linker_len)[0]} for p in range(4, 11)]
            c = alt.Chart(pd.DataFrame(res)).mark_line(point=True).encode(x='pH', y='Score').properties(height=200)
            st.altair_chart(c, use_container_width=True)
            if warnings:
                st.error(f"{len(warnings)} Critical Flags")
                for w in warnings: st.caption(f"🚨 {w}")
            for r in reasons: st.success(r)

        with tab_crowd:
            st.markdown("#### Lattice Simulator")
            spacing, coverage, status, crowd_warn = lab_brain.calculate_crowding_metrics(df_surface, target_density)
            c1, c2 = st.columns(2)
            c1.metric("Spacing", f"{spacing:.1f} Å")
            c2.metric("Coverage", f"{coverage:.0f}%")
            if "CRITICAL" in str(crowd_warn): st.error(status)
            elif "Risk" in str(crowd_warn): st.warning(status)
            else: st.success(status)
            for w in crowd_warn: st.caption(f"⚠️ {w}")

        with tab_actions:
            st.markdown("#### AI Architect")
            if not api_key: st.warning("OpenAI Key required.")
            elif not ai_suggestions: st.info("No mutations recommended.")
            else:
                for sugg in ai_suggestions:
                    with st.container(border=True):
                        st.markdown(f"**Suggestion:** `{sugg['mutation_code']}`")
                        st.caption(sugg['llm_reasoning'])
                        if st.button("Apply Mutation", key=sugg['mutation_code']):
                            parts = sugg['mutation_code'].split()
                            res_name, res_id_str = parts[0].split('-')
                            new_pdb = mutate_residue(pdb_str, f"{res_name}-{res_id_str}", parts[2])
                            if new_pdb:
                                st.session_state['mutated_pdb_str'] = new_pdb
                                st.rerun()

        with tab_calib:
            st.markdown("#### Lab-in-the-Loop")
            st.caption("Teach the model. Input your real lab results for this protein/surface pair.")
            
            c1, c2 = st.columns(2)
            c1.metric("Predicted Score", f"{score:.0f}")
            real_score = c2.number_input("Actual Lab Yield (%)", 0, 100, int(score))
            
            if st.button("➕ Add to Training Set"):
                datapoint = {
                    "pdb": pdb_input, "surface": s_choice,
                    "components": components, 
                    "actual": real_score
                }
                st.session_state['calibration_data'].append(datapoint)
                st.success(f"Added Data Point #{len(st.session_state['calibration_data'])}")
            
            if st.session_state['calibration_data']:
                st.markdown("---")
                st.caption(f"Training Set: {len(st.session_state['calibration_data'])} points")
                st.dataframe(pd.DataFrame(st.session_state['calibration_data']).drop(columns=['components']), hide_index=True)
                
                if st.button("Run Bayesian Refinement", type="primary"):
                    with st.spinner("Optimizing Physics Engine..."):
                        result_msg = lab_brain.optimize_weights(st.session_state['calibration_data'])
                        st.success(result_msg)
                        time.sleep(1)
                        st.rerun()
            else:
                st.info("Add at least 2 data points to calibrate.")
                
            st.markdown("---")
            st.caption("Current Physics Weights")
            st.json(lab_brain.weights)

        with tab_export:
            st.markdown("#### Manufacturing")
            audit_txt = generate_audit_report(pdb_input, s_choice, linker_label, ph, score, warnings)
            st.download_button("Compliance Report", audit_txt, "audit_report.txt", use_container_width=True)
            bot_script = generate_opentrons_protocol(pdb_input, linker_label, ph)
            st.download_button("Opentron Protocol code (.py)", bot_script, f"opentrons_{pdb_input}.py", use_container_width=True)