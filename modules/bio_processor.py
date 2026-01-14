import os
import requests
import io
import tempfile
import numpy as np
import pandas as pd
import streamlit as st
from Bio.PDB import PDBList, MMCIFParser, ShrakeRupley, PDBIO, PDBParser
from Bio.PDB.Polypeptide import is_aa
from modules.evolutionary_engine import EvolutionaryEngine

HYDROPHOBICITY_SCALE = {
    'ILE': 0.73, 'PHE': 0.61, 'VAL': 0.54, 'LEU': 0.53, 'TRP': 0.37,
    'MET': 0.26, 'ALA': 0.25, 'GLY': 0.16, 'CYS': 0.04, 'TYR': 0.02,
    'PRO': -0.07, 'THR': -0.18, 'SER': -0.26, 'HIS': -0.40, 'GLU': -0.62,
    'ASN': -0.64, 'GLN': -0.69, 'ASP': -0.72, 'LYS': -1.10, 'ARG': -1.76
}

AA3_TO_AA1 = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
    'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
    'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G'
}

def three_to_one(resname):
    """Safely converts 3-letter code to 1-letter code."""
    return AA3_TO_AA1.get(resname, 'X')

# --- HEURISTIC FALLBACK (Only used if API fails) ---
def _calculate_heuristic(residue, min_dist, sasa):
    PRIORS = {'CYS': 0.9, 'TRP': 0.9, 'HIS': 0.8, 'MET': 0.7}
    seq_score = PRIORS.get(residue.resname, 0.3)
    func_score = 1.0 if min_dist < 8.0 else 0.0
    struct_score = 1.0 if sasa < 15.0 else 0.2
    return (seq_score * 0.3) + (struct_score * 0.3) + (func_score * 0.4)


@st.cache_data(show_spinner=False)
def run_evolutionary_analysis(pdb_id, _structure):
    """
    Wrapper to run the EvolutionaryEngine with Caching.
    Passing _structure (underscore) prevents streamlit from hashing the whole object.
    We reconstruct the sequence string to hash.
    """
    # Extract Sequence from Structure
    query_seq = ""
    res_indices = [] # keep track of PDB IDs to map back
    
    for model in _structure:
        for chain in model:
            for residue in chain:
                if is_aa(residue, standard=True):
                    query_seq += three_to_one(residue.resname)
                    res_indices.append(residue.id[1])
            break # Process first chain only
        break # Process first model only
        
    engine = EvolutionaryEngine()
    cons_map_raw = engine.get_conservation_map(pdb_id, query_seq)
    
    if not cons_map_raw: return None
    
    # Remap 0-based index to PDB Residue ID
    final_map = {}
    for i, score in cons_map_raw.items():
        if i < len(res_indices):
            pdb_res_id = res_indices[i]
            final_map[pdb_res_id] = score
            
    return final_map

def _analyze_structure(structure, zinc_atoms, pdb_id):
    sr = ShrakeRupley()
    try: sr.compute(structure, level="R")
    except: pass 

    real_cons_map = run_evolutionary_analysis(pdb_id, structure)
    using_api = real_cons_map is not None

    data = []
    zn_vecs = [np.array(z.get_coord()) for z in zinc_atoms]
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if not is_aa(residue, standard=True): continue
                
                sasa = getattr(residue, "sasa", 0.0)
                coords = np.array([0.,0.,0.])
                if "CA" in residue: coords = residue["CA"].get_coord()
                
                min_dist = 999.0
                if zn_vecs and "CA" in residue:
                    dists = [np.linalg.norm(coords - z) for z in zn_vecs]
                    min_dist = min(dists)
                
                if using_api and residue.id[1] in real_cons_map:
                    cons_score = real_cons_map[residue.id[1]] / 100.0
                    method = "RCSB_Entropy"
                else:
                    cons_score = _calculate_heuristic(residue, min_dist, sasa)
                    method = "Struct_Heuristic"

                # Normalize 0-100
                cons_percent = min(100, max(0, cons_score * 100))
                
                h_raw = HYDROPHOBICITY_SCALE.get(residue.resname, 0.0)
                h_score = (h_raw + 1.76) / 2.49 * 100
                h_score = max(0, min(100, h_score))

                for atom in residue: atom.set_bfactor(cons_percent)

                data.append({
                    "Residue": residue.resname,
                    "Chain": chain.id,
                    "ID": residue.id[1],
                    "SASA": round(sasa, 1),
                    "Dist_to_Zinc": round(min_dist, 1),
                    "Charge": _get_charge(residue.resname),
                    "Conservation": round(cons_score, 2),
                    "Method": method,
                    "Hydrophobicity": round(h_score, 1),
                    "coords": coords 
                })

    df_surface = pd.DataFrame(data)
    
    io_w = PDBIO()
    io_w.set_structure(structure)
    stream = io.StringIO()
    io_w.save(stream)
    pdb_string = stream.getvalue()
    
    return df_surface, pdb_string

def _get_charge(resname):
    if resname in ['ARG', 'LYS']: return 1
    if resname in ['ASP', 'GLU']: return -1
    if resname == 'HIS': return 0.1 
    return 0

def _fetch_metadata(pdb_id):
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
    except: pass
    return {"title": "Unknown", "resolution": 9.99, "method": "Unknown"}

def map_metric_to_bfactor(pdb_string, df_surface, metric_col):
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("viz", stream)
    metric_map = dict(zip(df_surface['ID'], df_surface[metric_col]))
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] in metric_map:
                    val = metric_map[residue.id[1]]
                    for atom in residue: atom.set_bfactor(val)
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
    try: structure = parser.get_structure(pdb_id, cif_path)
    except: 
        fallback = os.path.join(temp_dir, f"{pdb_id}.cif")
        structure = parser.get_structure(pdb_id, fallback)

    zinc_atoms = [atom for atom in structure.get_atoms() if atom.element == 'ZN']
    
    df_surface, pdb_string = _analyze_structure(structure, zinc_atoms, pdb_id)
    zinc_coords = [z.get_coord().tolist() for z in zinc_atoms]
    
    return df_surface, pdb_string, zinc_coords, _fetch_metadata(pdb_id)

def process_pdb_string(pdb_string: str):
    parser = PDBParser(QUIET=True)
    stream = io.StringIO(pdb_string)
    structure = parser.get_structure("modified", stream)
    zinc_atoms = [atom for atom in structure.get_atoms() if atom.element == 'ZN']
    df_surface, _ = _analyze_structure(structure, zinc_atoms, "USER_UPLOAD")
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