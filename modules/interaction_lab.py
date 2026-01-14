import pandas as pd
import numpy as np
from modules.llm_factory import LLMFactory

class InteractionLab:
    """
    Reasoning Engine v9.0 - Aggressive Problem Detection.
    """
    def __init__(self):
        self.weights = {"electrostatic": 1.0, "hydrophobic": 1.0, "steric": 1.0}

    def calibrate(self, predicted_score, real_lab_score):
        delta = real_lab_score - predicted_score
        if delta < -10: 
            self.weights["steric"] *= 1.2
            self.weights["electrostatic"] *= 1.1
        elif delta > 10: 
            self.weights["hydrophobic"] *= 1.1
        return self.weights

    def calculate_interface_score(self, protein_features: pd.DataFrame, surface_props: dict, 
                                  rotation_matrix=None, ph: float = 7.4, 
                                  density: str = "Single", topology: str = "Flat",
                                  api_key: str = None):
        
        # Initialize Brain
        llm_brain = LLMFactory(api_key)
        
        score = 50 
        reasons = []
        warnings = []
        ai_suggestions = []
        
        if surface_props is None or protein_features.empty:
            return 0, ["Insufficient data."], [], []

        # 1. Apply Rotation
        coords = np.vstack(protein_features['coords'].values)
        if rotation_matrix is not None:
            com = np.mean(coords, axis=0)
            coords = (coords - com) @ rotation_matrix.T + com
            
        # 2. Physics & Logic Gates (Kill Zones)
        zinc_mask = protein_features['Dist_to_Zinc'] < 0.1 
        if not zinc_mask.all() and zinc_mask.any():
             zinc_indices = protein_features[zinc_mask].index
             min_zinc_z = np.min(coords[zinc_indices, 2])
             if min_zinc_z < -8.0:
                 return 10, [], [f"CRITICAL: Active Site Crash (Zinc Z={min_zinc_z:.1f})."], []

        z_vals = coords[:, 2]
        interface_mask = z_vals < -5.0
        interface_residues = protein_features[interface_mask].copy()
        
        if interface_residues.empty:
            return 20, ["Poor Orientation."], ["Optimize orientation."], []

        # 3. AI Analysis - EXPANDED LOGIC
        # Determine Surface Nature
        tpsa = surface_props.get('tpsa', 0)
        is_polar_surface = tpsa > 40
        is_positive_surface = "positive" in surface_props.get('charge_density', '').lower()

        for idx, res in interface_residues.iterrows():
            c = res['Charge']
            res_name = res['Residue']
            
            # pH logic
            if res_name == 'HIS' and ph > 6.0: c = 0
            if res_name == 'CYS' and ph > 8.0: c = -1
            if res_name in ['LYS', 'ARG'] and ph > 10.0: c = 0
            
            is_hydrophobic_res = res_name in ['LEU','ILE','VAL','PHE','TRP', 'MET', 'ALA']
            is_charged_res = c != 0
            
            problem = None
            
            # --- TRIGGER 1: Steric Risk (Bulky residues near surface) ---
            # If a large residue is very close to surface plane (Z < -9.0)
            # Z coords are in coords array, need to map back via index
            # Simplified: checking residue type + interface status
            if res_name in ['TRP', 'TYR', 'PHE', 'ARG'] and idx in interface_residues.index:
                 # Check actual Z depth
                 atom_z = coords[idx, 2]
                 if atom_z < -9.5: # Digging deep into surface
                     problem = "Steric Clash (Bulky Sidechain)"

            # --- TRIGGER 2: Hydrophobic Mismatch ---
            if not problem:
                if is_polar_surface and is_hydrophobic_res:
                    problem = "Solubility Mismatch (Hydrophobic on Polar)"
                elif not is_polar_surface and is_charged_res:
                    problem = "Desolvation Penalty (Charge on Hydrophobic)"

            # --- TRIGGER 3: Electrostatic Clash ---
            if not problem and is_positive_surface and c > 0:
                 problem = "Electrostatic Repulsion (Cation-Cation)"

            # Ask LLM (Limit to 2 suggestions to save time/cost)
            if problem and len(ai_suggestions) < 2:
                insight = llm_brain.generate_suggestion(res_name, res['ID'], problem, "Polar" if is_polar_surface else "Hydrophobic")
                # Deduplicate: Don't suggest same mutation twice
                if not any(s['mutation_code'] == insight['mutation_code'] for s in ai_suggestions):
                    ai_suggestions.append(insight)

        # 4. Standard Scoring
        local_charge = interface_residues['Charge'].sum()
        hydrophobic_count = interface_residues['Residue'].isin(['LEU','ILE','VAL','PHE','TRP']).sum()
        
        if is_polar_surface:
            if local_charge < -1: 
                 score += 30 * self.weights["electrostatic"]
                 reasons.append("Electrostatic Lock")
            elif local_charge > 1:
                 score -= 20 * self.weights["electrostatic"]
                 warnings.append("Electrostatic Repulsion")
        else:
            if hydrophobic_count > 5:
                score += 35 * self.weights["hydrophobic"]
                reasons.append("Hydrophobic Anchor")

        if density == "High (Monolayer)" and len(interface_residues) > 30:
            score -= 20 * self.weights["steric"]
            warnings.append("Crowding Penalty")

        return max(0, min(100, score)), reasons, warnings, ai_suggestions