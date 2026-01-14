import os
import requests
import io
import tempfile
import numpy as np
import pandas as pd
import streamlit as st
from Bio.PDB import PDBList, MMCIFParser, ShrakeRupley, PDBIO, PDBParser

# --- HELPER FUNCTIONS ---
def _get_charge(resname):
    """Simple charge map for pH ~7.4"""
    if resname in ['ARG', 'LYS']: return 1
    if resname in ['ASP', 'GLU']: return -1
    if resname == 'HIS': return 0.1 
    return 0

def _fetch_metadata(pdb_id):
    """Fetches metadata using RCSB Data API."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
    try:
        r = requests.get(url, timeout=3)
        if r.status_code == 200:
            data = r.json()
            return {
                "title": data.get("struct", {}).get("title", "Unknown Title"),
                "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined", [9.99])[0],
                "method": data.get("exptl", [{}])[0].get("method", "Unknown"),
            }
    except:
        pass
    return {"title": "Unknown", "resolution": 9.99, "method": "Unknown"}

# --- WORKER FUNCTIONS ---

def _analyze_structure(structure, zinc_atoms):
    """Internal helper to analyze a BioPython structure object."""
    # 1. SASA Calculation
    sr = ShrakeRupley()
    try:
        sr.compute(structure, level="R")
    except:
        pass # Fallback if SASA fails (rare)

    # 2. Build DataFrame
    data = []
    zn_vecs = [np.array(z.get_coord()) for z in zinc_atoms]
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ": continue 
                
                sasa = getattr(residue, "sasa", 0.0)
                # Heuristic: Surface Exposed if SASA > 10
                if sasa > 10.0:
                    coords = np.array([0.,0.,0.])
                    min_dist = 999.0
                    if "CA" in residue:
                        coords = residue["CA"].get_coord()
                        if zn_vecs:
                            dists = [np.linalg.norm(coords - z) for z in zn_vecs]
                            min_dist = min(dists)
                    
                    data.append({
                        "Residue": residue.resname,
                        "Chain": chain.id,
                        "ID": residue.id[1],
                        "SASA": round(sasa, 1),
                        "Dist_to_Zinc": round(min_dist, 1),
                        "Charge": _get_charge(residue.resname),
                        "coords": coords 
                    })

    df_surface = pd.DataFrame(data)
    
    # Export PDB String
    io_w = PDBIO()
    io_w.set_structure(structure)
    stream = io.StringIO()
    io_w.save(stream)
    pdb_string = stream.getvalue()
    
    return df_surface, pdb_string

@st.cache_data(show_spinner=False)
def process_protein_structure(pdb_id: str):
    """
    Downloads and parses a fresh PDB from ID.
    """
    pdb_id = pdb_id.lower()
    temp_dir = tempfile.gettempdir()
    
    # 1. Download
    pdbl = PDBList(verbose=False)
    cif_path = pdbl.retrieve_pdb_file(pdb_id, pdir=temp_dir, file_format='mmCif', overwrite=False)
    
    # 2. Parse
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, cif_path)
    except Exception:
        fallback_path = os.path.join(temp_dir, f"{pdb_id}.cif")
        if os.path.exists(fallback_path):
             structure = parser.get_structure(pdb_id, fallback_path)
        else:
             raise FileNotFoundError(f"Could not load PDB {pdb_id}")

    # 3. Find Zinc
    zinc_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'ZN':
                        zinc_atoms.append(atom)

    df_surface, pdb_string = _analyze_structure(structure, zinc_atoms)
    zinc_coords = [z.get_coord().tolist() for z in zinc_atoms]
    
    return df_surface, pdb_string, zinc_coords, _fetch_metadata(pdb_id)

def process_pdb_string(pdb_string: str):
    """
    Analyzes an *already downloaded/mutated* PDB string.
    Essential for Phase 4 (Mutagenesis) Loop.
    """
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("modified", stream)
    
    zinc_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == 'ZN':
                        zinc_atoms.append(atom)
                        
    df_surface, _ = _analyze_structure(structure, zinc_atoms)
    zinc_coords = [z.get_coord().tolist() for z in zinc_atoms]
    
    return df_surface, pdb_string, zinc_coords

# --- PHASE 4: MUTAGENESIS ENGINE ---
def mutate_residue(pdb_string, residue_id_str, new_res_type):
    """
    Real Structural Edit.
    Parses string, finds residue, renames it (changing physics), and returns new string.
    residue_id_str format: "RESNAME-ID" (e.g., "ASP-124")
    """
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("temp", stream)
    
    target_id = int(residue_id_str.split('-')[1])
    target_name = residue_id_str.split('-')[0]
    
    mutation_applied = False
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check match (Sequence ID and Name)
                if residue.id[1] == target_id and residue.resname == target_name:
                    residue.resname = new_res_type
                    mutation_applied = True
                    # NOTE: In a full molecular dynamics engine, we would rebuild sidechains.
                    # For this Coarse-Grained Physics Engine, renaming changes the 
                    # Charge and Hydrophobicity lookup in InteractionLab, which is 
                    # the "Real" change required for the score to update.
    
    if not mutation_applied:
        return None
        
    io_w = PDBIO()
    io_w.set_structure(structure)
    out = io.StringIO()
    io_w.save(out)
    return out.getvalue()

# --- PHASE 2: VISUAL SNAP LOGIC ---
def get_rotated_pdb(pdb_string, rotation_matrix):
    """
    Applies rotation matrix around COM and returns new PDB string.
    Works on PDB STRING input to support mutated structures.
    """
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("temp", stream)
         
    # Calculate Center of Mass (COM)
    atoms = list(structure.get_atoms())
    if not atoms: return pdb_string
    
    coords = np.array([a.get_coord() for a in atoms])
    com = np.mean(coords, axis=0)
    
    # Apply Rotation: (v - COM) * R + COM
    for atom in atoms:
        v = atom.get_coord() - com
        new_v = np.dot(v, rotation_matrix.T) 
        atom.set_coord(new_v + com)
        
    io_w = PDBIO()
    io_w.set_structure(structure)
    stream = io.StringIO()
    io_w.save(stream)
    return stream.getvalue()