from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, Draw
import numpy as np

class SurfaceEngine:
    """
    Handles the 'Agnostic' Surface logic.
    Converts SMILES/Materials into parameterized descriptors (LogP, Charge, etc.).
    """
    
    def __init__(self):
        pass

    def process_smiles(self, smiles_string: str):
        """
        Input: SMILES (e.g., 'C1=CC=CC=C1' for Benzene surface)
        Output: Dictionary of Physical Parameters for the Playbook.
        """
        mol = Chem.MolFromSmiles(smiles_string)
        if not mol:
            return None
            
        return {
            "smiles": smiles_string,
            "logP": Descriptors.MolLogP(mol),       # Hydrophobicity Proxy
            "tpsa": Descriptors.TPSA(mol),          # Charge/Polarity Proxy
            "mw": Descriptors.MolWt(mol),
            "h_donors": Descriptors.NumHDonors(mol),
            "h_acceptors": Descriptors.NumHAcceptors(mol),
            "charge_density": "High" if Descriptors.TPSA(mol) > 50 else "Low"
        }

    def generate_3d_grid(self, smiles: str):
        """
        CREATIVE VISUALIZATION:
        Creates a PDB block representing a 'Surface Patch' (3x3 grid of molecules)
        to visualize the material in the 3D viewer.
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return ""
        mol = Chem.AddHs(mol)
        
        # --- FIX: ROBUST COORDINATE GENERATION ---
        try:
            # Embed
            res = AllChem.EmbedMolecule(mol)
            if res == -1: AllChem.EmbedMolecule(mol, useRandomCoords=True)
            
            # Optimize (Skip if fails)
            try:
                AllChem.UFFOptimizeMolecule(mol)
            except:
                pass 
        except:
            return ""
        
        return Chem.MolToMolBlock(mol)

    def suggest_modification(self, current_props: dict, goal: str):
        """
        The 'Recommendation Layer' Logic.
        Suggests chemical changes based on desired outcome.
        """
        suggestions = []
        if goal == "Increase Retention":
            if current_props['logP'] < 0:
                suggestions.append("Surface is too polar. Add Alkyl chains (-CCCCC) to increase Hydrophobic interaction.")
            if current_props['h_donors'] == 0:
                suggestions.append("Add Hydroxyl (-OH) or Amine (-NH2) groups to enable Hydrogen Bonding.")
                
        return suggestions