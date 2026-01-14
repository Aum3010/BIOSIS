import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull
from scipy.optimize import minimize
from modules.llm_factory import LLMFactory

class InteractionLab:
    
    def __init__(self):
        # Default Weights (The "Prior")
        # These will be updated by the optimizer
        self.weights = {
            "electrostatic": 1.0, 
            "hydrophobic": 1.0, 
            "steric": 1.0
        }

    def _get_dynamic_charge(self, res_name, ph):
        if res_name in ['ASP', 'GLU']: return -1 if ph > 4.0 else 0
        if res_name == 'CYS': return -1 if ph > 8.3 else 0
        if res_name == 'TYR': return -1 if ph > 10.0 else 0
        if res_name == 'HIS': return 1 if ph < 6.0 else 0
        if res_name == 'LYS': return 1 if ph < 10.5 else 0
        if res_name == 'ARG': return 1 
        return 0

    def calculate_crowding_metrics(self, protein_features, density_pmol_cm2):
        if protein_features.empty: return 0, 0, "No Data", []
        
        coords = np.vstack(protein_features['coords'].values)
        xy_points = coords[:, :2]
        
        try:
            hull = ConvexHull(xy_points)
            footprint_area = hull.volume
        except: 
            footprint_area = 200.0
            
        radius = np.sqrt(footprint_area / np.pi)
        diameter = radius * 2

        if density_pmol_cm2 <= 0.01: 
            return 999.0, 0.0, "Free Diffusion", []

        molecules_per_cm2 = density_pmol_cm2 * 1e-12 * 6.022e23
        area_per_molecule_A2 = 1e16 / molecules_per_cm2
        spacing = np.sqrt(area_per_molecule_A2 / 0.866)
        coverage = min(100, (footprint_area / area_per_molecule_A2) * 100)
        
        status = "🟢 Free Diffusion"
        warnings = []
        
        if spacing < diameter:
            status = "🔴 Steric Clash"
            warnings.append(f"CRITICAL: Proteins overlap! (Spacing {spacing:.1f}Å < Width {diameter:.1f}Å)")
        elif spacing < (diameter + 10):
            status = "🟡 Jamming Limit"
            warnings.append("Crowding Risk: Diffusion limited.")
            
        return spacing, coverage, status, warnings

    def calculate_interface_score(self, protein_features: pd.DataFrame, surface_props: dict, 
                                  rotation_matrix=None, ph: float = 7.4, 
                                  density: str = "Single", topology: str = "Flat",
                                  api_key: str = None,
                                  linker_length: float = 0.0,
                                  linker_type: str = "flexible"):
        
        llm_brain = LLMFactory(api_key)
        
        score = 50 
        reasons = []
        warnings = []
        ai_suggestions = []
        
        components = {"electrostatic": 0.0, "hydrophobic": 0.0, "steric": 0.0}
        
        if surface_props is None or protein_features.empty:
            return 0, ["Insufficient data."], [], [], components

        coords = np.vstack(protein_features['coords'].values)
        if rotation_matrix is not None:
            com = np.mean(coords, axis=0)
            coords = (coords - com) @ rotation_matrix.T + com
            
        zinc_mask = protein_features['Dist_to_Zinc'] < 0.1 
        if not zinc_mask.all() and zinc_mask.any():
             zinc_indices = protein_features[zinc_mask].index
             min_zinc_z = np.min(coords[zinc_indices, 2])
             if min_zinc_z < -5.0:
                 if linker_length < 10.0:
                     return 10, [], [f"CRITICAL: Active Site Blocked."], [], components
                 else:
                     penalty = 20
                     score -= penalty * self.weights['steric']
                     components['steric'] -= penalty
                     warnings.append("Orientation Risk: Active site faces surface.")

        if linker_length > 0:
            bonus = min(25, linker_length * 1.25)
            # Linker acts as a steric buffer
            score += bonus * self.weights['steric']
            components['steric'] += bonus
            reasons.append(f"Linker (+{int(bonus)}): Reduced Denaturation Risk")
            
            wobble = 0
            if linker_type == "flexible" and linker_length > 15:
                wobble = ((linker_length - 15) ** 1.5) * 0.5
            elif linker_type == "rigid" and linker_length > 30:
                wobble = 10
            
            if wobble > 0:
                score -= wobble * self.weights['steric']
                components['steric'] -= wobble
                warnings.append(f"Entropy Penalty: Linker instability.")

        bottom_threshold = -5.0 - linker_length
        z_vals = coords[:, 2]
        interface_mask = z_vals < (bottom_threshold + 5.0) 
        interface_residues = protein_features[interface_mask].copy()

        surf_logp = surface_props.get('logP', 0)
        is_surf_hydrophobic = surf_logp > 1.5
        
        if not interface_residues.empty:
            avg_int_hydro = interface_residues['Hydrophobicity'].mean()
            if pd.isna(avg_int_hydro): avg_int_hydro = 0
            
            if is_surf_hydrophobic and avg_int_hydro < 40:
                penalty = 25
                if linker_length > 10: penalty = 5 
                
                score -= penalty * self.weights['hydrophobic']
                components['hydrophobic'] -= penalty
                warnings.append(f"Denaturation Risk: Polar surface on Hydrophobic material.")
            elif is_surf_hydrophobic and avg_int_hydro > 60:
                bonus = 15
                score += bonus * self.weights['hydrophobic']
                components['hydrophobic'] += bonus
                reasons.append("Hydrophobic Anchor.")

        if linker_length == 0:
            if interface_residues.empty:
                 return 20, ["Poor Orientation."], ["Optimize orientation."], [], components
            
            surface_type = "neutral"
            if "positive" in surface_props.get('charge_density', '').lower(): surface_type = "positive"
            
            net_interface_charge = 0
            for idx, res in interface_residues.iterrows():
                c = self._get_dynamic_charge(res['Residue'], ph)
                net_interface_charge += c
                
                # Single residue logic
                problem = None
                if surface_type == "positive" and c > 0: problem = "Electrostatic Repulsion"
                
                if problem:
                    is_conserved = res.get('Conservation', 0.5) > 0.8
                    if is_conserved:
                        score -= 5 * self.weights['electrostatic']
                        components['electrostatic'] -= 5
                        warnings.append(f"Conserved Residue conflict.")
                    elif len(ai_suggestions) < 2:
                        insight = llm_brain.generate_suggestion(res['Residue'], res['ID'], problem, surface_type)
                        if not any(s['mutation_code'] == insight['mutation_code'] for s in ai_suggestions):
                            ai_suggestions.append(insight)
                            
            is_surface_polar = (surface_props.get('tpsa', 0) > 40)
            if is_surface_polar:
                val = 30
                if net_interface_charge < -2: 
                     score += val * self.weights["electrostatic"]
                     components['electrostatic'] += val
                     reasons.append(f"Electrostatic Lock")
                elif net_interface_charge > 2:
                     score -= val * self.weights["electrostatic"]
                     components['electrostatic'] -= val
                     warnings.append(f"Electrostatic Repulsion")

        return max(0, min(100, score)), reasons, warnings, ai_suggestions, components

    def optimize_weights(self, training_data):
        """
        Bayesian-style weight refinement using Least Squares Minimization.
        training_data: List of dicts [{'components': {...}, 'actual': 85}, ...]
        """
        if len(training_data) < 2: return "Need at least 2 data points."

        def loss_function(weights_array):
            w_e, w_h, w_s = weights_array
            error_sum = 0
            
            for point in training_data:
                comps = point['components']
                # Reconstruct score: Base(50) + Weighted Components
                pred = 50 + (comps['electrostatic'] * w_e) + \
                            (comps['hydrophobic'] * w_h) + \
                            (comps['steric'] * w_s)
                
                # Clamp to 0-100 logic
                pred = max(0, min(100, pred))
                error_sum += (pred - point['actual']) ** 2
                
            return error_sum

        x0 = [self.weights['electrostatic'], self.weights['hydrophobic'], self.weights['steric']]
        
        # Bounds (Prevent negative weights or infinite scaling)
        # We allow weights to vary from 0.1x to 5.0x
        bnds = ((0.1, 5.0), (0.1, 5.0), (0.1, 5.0))
        
        res = minimize(loss_function, x0, bounds=bnds, method='L-BFGS-B')
        
        if res.success:
            self.weights['electrostatic'] = round(res.x[0], 2)
            self.weights['hydrophobic'] = round(res.x[1], 2)
            self.weights['steric'] = round(res.x[2], 2)
            return f"Converged! New Weights: E={self.weights['electrostatic']}, H={self.weights['hydrophobic']}, S={self.weights['steric']}"
        else:
            return "Optimization failed. Data might be contradictory."