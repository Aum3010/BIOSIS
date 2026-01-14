import streamlit as st
import pandas as pd
import time
import hashlib
import json
import altair as alt
from modules.bio_processor import process_protein_structure, process_pdb_string, mutate_residue, get_rotated_pdb
from modules.surface_engine import SurfaceEngine
from modules.interaction_lab import InteractionLab
from modules.literature_scout import LiteratureScout
from modules.visualizer import SurfaceVisualizer
from modules.optimizer import GeneticOptimizer
from modules.generator import MaterialGenerator
from stmol import showmol
import py3Dmol

# --- CONFIG ---
st.set_page_config(page_title="BIOSIS Architect", layout="wide", page_icon="🧬")

@st.cache_resource
def get_static_engines():
    return SurfaceEngine(), InteractionLab(), LiteratureScout(), SurfaceVisualizer(), MaterialGenerator()

chem_engine, lab_brain, scout, vis_engine, generator = get_static_engines()

if 'current_pdb' not in st.session_state: st.session_state['current_pdb'] = "1V9E"
if 'mutated_pdb_str' not in st.session_state: st.session_state['mutated_pdb_str'] = None
if 'generated_surface' not in st.session_state: st.session_state['generated_surface'] = None

# --- SIDEBAR: SYSTEM DEFINITION ---
with st.sidebar:
    st.header("1. Configuration")
    api_key = st.text_input("OpenAI API Key", type="password", help="Required for AI Suggestions")
    
    st.subheader("Protein")
    pdb_input = st.text_input("PDB ID", st.session_state['current_pdb']).upper()
    if pdb_input != st.session_state['current_pdb']:
        st.session_state['current_pdb'] = pdb_input
        st.session_state['mutated_pdb_str'] = None
        st.rerun()

    if st.session_state['mutated_pdb_str']:
        st.warning("🧪 Variant Active")
        if st.button("Reset to Wild Type"):
            st.session_state['mutated_pdb_str'] = None
            st.rerun()

    st.subheader("Surface")
    s_choice = st.selectbox("Material", ["Graphene", "Silica", "Gold", "Custom"])
    base_smiles = {"Graphene": "c1ccccc1", "Silica": "O=[Si]=O", "Gold": "[Au]"}.get(s_choice, "c1ccccc1")
    if s_choice == "Custom": base_smiles = st.text_input("SMILES", "NCC")

    st.subheader("Reactor")
    ph = st.slider("pH", 4.0, 10.0, 7.4)
    density = st.select_slider("Crowding", ["Single", "High (Monolayer)"])
    topology = st.selectbox("Topology", ["Flat", "Nanopore", "Nanopillar"])
    
    st.markdown("---")
    run_sim = st.button("🚀 Run Simulation", type="primary", use_container_width=True)

# --- MAIN WORKBENCH ---
st.title("BIOSIS Architect v9.0")

if not run_sim and not st.session_state.get('last_sim'):
    st.info("👋 **Ready.** Configure your system in the sidebar and click **Run Simulation**.")
else:
    st.session_state['last_sim'] = True
    
    with st.spinner("Processing Digital Twin..."):
        # 1. Load Data
        if st.session_state['mutated_pdb_str']:
            df_surface, pdb_str, zinc_coords = process_pdb_string(st.session_state['mutated_pdb_str'])
            meta = {"title": "Engineered Variant", "method": "In-Silico"}
        else:
            df_surface, pdb_str, zinc_coords, meta = process_protein_structure(pdb_input)

        surface_props = chem_engine.process_smiles(base_smiles)
        
        # 2. Optimize
        optimizer = GeneticOptimizer(df_surface, zinc_coords)
        t_type = "positive" if surface_props['tpsa'] > 40 else "hydrophobic"
        best_rot, _, _ = optimizer.run(surface_type=t_type)
        pdb_rotated = get_rotated_pdb(pdb_str, best_rot)
        
        # 3. Score
        score, reasons, warnings, ai_suggestions = lab_brain.calculate_interface_score(
            df_surface, surface_props, best_rot, ph, density, topology, api_key
        )

    # --- ZONE 1: SCOREBOARD (Top) ---
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("BIOSIS Score", f"{score:.0f}", help="0-100 Success Probability")
    c2.metric("Surface Charge", surface_props['charge_density'])
    c3.metric("Orientation", "Optimized")
    c4.metric("Risk Assessment", "Low" if score > 70 else "High", delta_color="inverse")

    st.markdown("---")

    # --- ZONE 2 & 3: SPLIT VIEW ---
    col_left, col_right = st.columns([1.8, 1.2])

    # LEFT: VISUAL TWIN
    with col_left:
        st.subheader("3D Reality Simulator")
        view = py3Dmol.view(height=600)
        view.addModel(pdb_rotated, 'pdb')
        view.setStyle({'cartoon': {'color': 'white'}})
        # Red Zinc
        for z in zinc_coords:
             view.addSphere({'center':{'x':z[0],'y':z[1],'z':z[2]}, 'radius': 2.0, 'color': 'red'})
        # Surface
        surface_pdb = vis_engine.generate_topology(base_smiles, shape=topology, density=density)
        if surface_pdb:
            view.addModel(surface_pdb, 'pdb')
            view.setStyle({'model': -1}, {'stick': {}})
        view.zoomTo()
        showmol(view, height=600)

    # RIGHT: INTELLIGENCE DECK
    with col_right:
        # We use Tabs to organize the "Messy" parts logically
        tab_analysis, tab_actions, tab_export = st.tabs(["🧠 Intelligence", "🛠️ Actions", "📂 Export"])
        
        with tab_analysis:
            # 1. Literature Feed (MOVED HERE)
            scout.render_feed(meta.get('title', pdb_input), s_choice)
            
            st.markdown("#### Stability Profile")
            res = [{"pH": p, "Score": lab_brain.calculate_interface_score(df_surface, surface_props, best_rot, float(p), density, topology, api_key)[0]} for p in range(4, 11)]
            c = alt.Chart(pd.DataFrame(res)).mark_line(point=True).encode(x='pH', y='Score').properties(height=200)
            st.altair_chart(c, use_container_width=True)

            if warnings:
                st.error(f"{len(warnings)} Critical Flags Detected")
                for w in warnings: st.caption(f"🚨 {w}")

        with tab_actions:
            st.markdown("#### AI Architect")
            if not api_key:
                st.warning("Enter OpenAI API Key in Sidebar to enable AI.")
            elif not ai_suggestions:
                st.info("No mutations recommended by AI.")
            else:
                for sugg in ai_suggestions:
                    with st.container(border=True):
                        st.markdown(f"**Suggestion:** `{sugg['mutation_code']}`")
                        st.caption(sugg['llm_reasoning'])
                        if st.button("Apply Mutation", key=sugg['mutation_code']):
                            parts = sugg['mutation_code'].split()
                            res_name, res_id = parts[0].split('-')
                            new_pdb = mutate_residue(pdb_str, f"{res_name}-{res_id}", parts[2])
                            if new_pdb:
                                st.session_state['mutated_pdb_str'] = new_pdb
                                st.rerun()

        with tab_export:
            st.markdown("#### Manufacturing")
            st.download_button("📄 Download Compliance Report", f"Hash: {hashlib.sha256(str(score).encode()).hexdigest()}", "audit.txt", use_container_width=True)
            st.download_button("🤖 Download Robot Protocol", "import opentrons...", "bot.py", use_container_width=True)
            
            st.markdown("#### Calibration")
            real = st.number_input("Lab Result (0-100)", 0, 100, 0)
            if st.button("Update Physics Model"):
                lab_brain.calibrate(score, real)
                st.success("Model Weights Updated")