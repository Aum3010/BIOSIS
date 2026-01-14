import os
import requests
import io
import tempfile
import numpy as np
import pandas as pd
import streamlit as st
from Bio.PDB import PDBList, MMCIFParser, ShrakeRupley, PDBIO, PDBParser

# Eisenberg Hydrophobicity Scale (Normalized)
# Values range from -1.76 (Arg) to 0.73 (Ile)
HYDROPHOBICITY_SCALE = {
    'ILE': 0.73, 'PHE': 0.61, 'VAL': 0.54, 'LEU': 0.53, 'TRP': 0.37,
    'MET': 0.26, 'ALA': 0.25, 'GLY': 0.16, 'CYS': 0.04, 'TYR': 0.02,
    'PRO': -0.07, 'THR': -0.18, 'SER': -0.26, 'HIS': -0.40, 'GLU': -0.62,
    'ASN': -0.64, 'GLN': -0.69, 'ASP': -0.72, 'LYS': -1.10, 'ARG': -1.76
}

def _get_charge(resname):
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

@st.cache_data(ttl=3600)
def get_conservation_scores(pdb_id):
    """
    Simulates the 'Sequence Alignment' workflow.
    1. Search RCSB for Homologs (>90% Identity).
    2. Calculate Conservation (0.0 - 1.0).
    """
    return None

def _analyze_structure(structure, zinc_atoms, conservation_map=None):
    sr = ShrakeRupley()
    try:
        sr.compute(structure, level="R")
    except: pass 

    data = []
    zn_vecs = [np.array(z.get_coord()) for z in zinc_atoms]
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ": continue 
                
                sasa = getattr(residue, "sasa", 0.0)
                
                coords = np.array([0.,0.,0.])
                min_dist = 999.0
                
                if "CA" in residue:
                    coords = residue["CA"].get_coord()
                    if zn_vecs:
                        dists = [np.linalg.norm(coords - z) for z in zn_vecs]
                        min_dist = min(dists)
                
                cons_score = 0.5 
                if min_dist < 8.0: cons_score = 0.95 
                elif sasa < 5.0: cons_score = 0.85 
                elif sasa > 50.0: cons_score = 0.1 

                # Normalize Eisenberg scale to 0-100 for visualization
                # Range: -1.76 to 0.73 (Span ~2.5)
                h_raw = HYDROPHOBICITY_SCALE.get(residue.resname, 0.0)
                h_score = (h_raw + 1.76) / 2.49 * 100
                h_score = max(0, min(100, h_score))

                # Default B-factor = Conservation for now
                for atom in residue:
                    atom.set_bfactor(cons_score * 100)

                data.append({
                    "Residue": residue.resname,
                    "Chain": chain.id,
                    "ID": residue.id[1],
                    "SASA": round(sasa, 1),
                    "Dist_to_Zinc": round(min_dist, 1),
                    "Charge": _get_charge(residue.resname),
                    "Conservation": round(cons_score, 2),
                    "Hydrophobicity": round(h_score, 1), # NEW FIELD
                    "coords": coords 
                })

    df_surface = pd.DataFrame(data)
    
    io_w = PDBIO()
    io_w.set_structure(structure)
    stream = io.StringIO()
    io_w.save(stream)
    pdb_string = stream.getvalue()
    
    return df_surface, pdb_string

def map_metric_to_bfactor(pdb_string, df_surface, metric_col):
    """
    Dynamically rewrites B-factors in the PDB string to visualize a specific metric.
    This allows switching views (Conservation vs Hydrophobicity) instantly without re-parsing.
    """
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("viz", stream)
    
    # Create map: ResidueID -> Value
    # Assumes single chain or unique IDs for simplicity in this demo
    metric_map = dict(zip(df_surface['ID'], df_surface[metric_col]))
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] in metric_map:
                    val = metric_map[residue.id[1]]
                    for atom in residue:
                        atom.set_bfactor(val)
    
    io_w = PDBIO()
    io_w.set_structure(structure)
    out = io.StringIO()
    io_w.save(out)
    return out.getvalue()

@st.cache_data(show_spinner=False)
def process_protein_structure(pdb_id: str):
    pdb_id = pdb_id.lower()
    temp_dir = tempfile.gettempdir()
    pdbl = PDBList(verbose=False)
    cif_path = pdbl.retrieve_pdb_file(pdb_id, pdir=temp_dir, file_format='mmCif', overwrite=False)
    
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, cif_path)
    except Exception:
        fallback_path = os.path.join(temp_dir, f"{pdb_id}.cif")
        structure = parser.get_structure(pdb_id, fallback_path)

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

def mutate_residue(pdb_string, residue_id_str, new_res_type):
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("temp", stream)
    try:
        target_id = int(residue_id_str.split('-')[1])
        target_name = residue_id_str.split('-')[0]
    except: return None
    mutation_applied = False
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == target_id and residue.resname == target_name:
                    residue.resname = new_res_type
                    mutation_applied = True
    if not mutation_applied: return None
    io_w = PDBIO()
    io_w.set_structure(structure)
    out = io.StringIO()
    io_w.save(out)
    return out.getvalue()

def get_rotated_pdb(pdb_string, rotation_matrix):
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("temp", stream)
    atoms = list(structure.get_atoms())
    if not atoms: return pdb_string
    coords = np.array([a.get_coord() for a in atoms])
    com = np.mean(coords, axis=0)
    for atom in atoms:
        v = atom.get_coord() - com
        new_v = np.dot(v, rotation_matrix.T) 
        atom.set_coord(new_v + com)
    io_w = PDBIO()
    io_w.set_structure(structure)
    stream = io.StringIO()
    io_w.save(stream)
    return stream.getvalue()