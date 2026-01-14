import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation as R

class GeneticOptimizer:
    """
    Evolutionary Algorithm to optimize Protein Orientation.
    Variables: Euler Angles (x, y, z)
    Fitness: Maximizing Surface Contact of 'Anchors' while keeping Zinc 'Safe'.
    """

    def __init__(self, protein_df: pd.DataFrame, zinc_coords_input):
        # 1. Extract Anchor Coordinates
        anchors_df = protein_df[protein_df['SASA'] > 10.0]
        if not anchors_df.empty and 'coords' in anchors_df.columns:
            self.anchor_coords = np.vstack(anchors_df['coords'].values)
            self.anchor_charges = anchors_df['Charge'].values
        else:
            self.anchor_coords = np.array([])
            self.anchor_charges = np.array([])

        # 2. Extract Zinc Coordinates (Handle input as list of lists OR list of objects)
        if zinc_coords_input and len(zinc_coords_input) > 0:
            # Check if input is list of lists (from cached JSON) or Atom objects
            if hasattr(zinc_coords_input[0], 'get_coord'):
                 self.zinc_coords = np.array([z.get_coord() for z in zinc_coords_input])
            else:
                 self.zinc_coords = np.array(zinc_coords_input)
        else:
            self.zinc_coords = np.array([])
        
        # 3. Calculate Center of Mass (COM)
        all_points = []
        if len(self.anchor_coords) > 0: all_points.append(self.anchor_coords)
        if len(self.zinc_coords) > 0: all_points.append(self.zinc_coords)
        
        if len(all_points) > 0:
            self.com = np.mean(np.vstack(all_points), axis=0)
        else:
            self.com = np.array([0.0, 0.0, 0.0])

    def fitness(self, rotation_matrix, target_charge_type):
        """
        Fitness Function = (Electrostatic Score) - (Steric Penalty)
        """
        score = 0
        
        # VIRTUAL SURFACE PLANE is at Z = -10 (relative to protein COM being at 0)
        SURFACE_Z = -10.0
        INTERACTION_ZONE = -5.0 # Region where binding happens
        
        # --- 1. Apply Rotation ---
        # Shift to origin -> Rotate -> Shift back
        # Note: We check positions relative to COM, assuming COM is at (0,0,0) for simulation
        
        # Rotate Zinc
        if len(self.zinc_coords) > 0:
            z_centered = self.zinc_coords - self.com
            z_rotated = z_centered @ rotation_matrix.T
            
            # CRITICAL CONSTRAINT: Zinc must point AWAY from surface
            # If Zinc is the lowest point, it hits the surface -> PENALTY
            min_z = np.min(z_rotated[:, 2])
            if min_z < -5.0: # If Zinc is pointing down
                return -500 # Heavy Penalty (Kill this solution)

        # Rotate Anchors
        if len(self.anchor_coords) > 0:
            a_centered = self.anchor_coords - self.com
            a_rotated = a_centered @ rotation_matrix.T
            
            z_vals = a_rotated[:, 2]
            
            # Find residues touching the virtual surface (bottom 5 Angstroms)
            interacting_indices = np.where(z_vals < -5.0)[0]
            
            if len(interacting_indices) == 0:
                return -50 # No contact with surface -> Bad orientation
            
            # Calculate Electrostatic Score for touching residues
            charges = self.anchor_charges[interacting_indices]
            
            if target_charge_type == "positive": # Surface is +
                score += np.sum(charges == -1) * 10  # Reward Negatives (Asp, Glu)
                score -= np.sum(charges == 1) * 10   # Punish Positives
            elif target_charge_type == "negative": # Surface is -
                score += np.sum(charges == 1) * 10   # Reward Positives (Lys, Arg)
                score -= np.sum(charges == -1) * 10
            elif target_charge_type == "hydrophobic":
                score += len(interacting_indices) * 2 # Reward just contact (VdW)
                
        return score

    def run(self, surface_type="neutral", generations=15, pop_size=30):
        """
        Runs the GA. Returns best rotation matrix and history.
        """
        best_score = -9999
        best_matrix = np.eye(3)
        history = []
        
        # Determine target
        target = "hydrophobic"
        if "polar" in surface_type or "positive" in surface_type: target = "positive"
        if "negative" in surface_type: target = "negative"

        # Initialize Population (Random Euler Angles)
        population = np.random.uniform(0, 360, (pop_size, 3))
        
        for g in range(generations):
            gen_scores = []
            
            # Evaluate
            for i in range(pop_size):
                r = R.from_euler('xyz', population[i], degrees=True)
                matrix = r.as_matrix()
                fit = self.fitness(matrix, target)
                gen_scores.append(fit)
                
                if fit > best_score:
                    best_score = fit
                    best_matrix = matrix
            
            history.append(max(gen_scores))
            
            # Selection & Mutation (Simple Tournament)
            # Take top 50% and mutate slightly
            sorted_idx = np.argsort(gen_scores)[::-1]
            top_half = population[sorted_idx[:pop_size//2]]
            
            # Create next gen: Top half + Mutated copies
            mutation_noise = np.random.normal(0, 15, top_half.shape) # +/- 15 degrees mutation
            next_gen = np.vstack([top_half, top_half + mutation_noise])
            
            population = next_gen[:pop_size] # Keep size constant

        return best_matrix, best_score, history