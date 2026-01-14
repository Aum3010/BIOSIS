from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import math

class SurfaceVisualizer:
    def generate_topology(self, smiles: str, shape: str = "Flat", density: str = "Single") -> str:
        """Generates PDB block for Surface Topology (Flat, Pore, Pillar)."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return ""
        mol = Chem.AddHs(mol)
        
        # --- FIX: ROBUST COORDINATE GENERATION ---
        try:
            # 1. Generate initial 3D coords (Random/Distance geometry)
            res = AllChem.EmbedMolecule(mol)
            if res == -1: # Try random coords if standard embed fails
                AllChem.EmbedMolecule(mol, useRandomCoords=True)
            
            # 2. Try to optimize, but ignore if it fails (e.g., for Gold [Au])
            try:
                AllChem.UFFOptimizeMolecule(mol)
            except Exception:
                pass # UFF failed (likely no params for atom), keeping raw coords
                
        except Exception:
            return ""

        conf = mol.GetConformer()
        atoms = mol.GetAtoms()
        pdb_lines = []
        atom_serial = 1
        res_serial = 1
        
        # Grid settings (Phase 7)
        grid_dim = 12 if density == "High (Monolayer)" else 6
        spacing = 5.0 
        
        for row in range(grid_dim):
            for col in range(grid_dim):
                # Phase 8: Geometry Logic
                if shape == "Flat":
                    x_base = (row * spacing) - ((grid_dim * spacing) / 2)
                    y_base = (col * spacing) - ((grid_dim * spacing) / 2)
                    z_base = -15.0
                elif shape == "Nanopore":
                    radius = 25.0
                    theta = (col / grid_dim) * 2 * math.pi
                    height = (row * spacing) - ((grid_dim * spacing) / 2)
                    x_base = radius * math.cos(theta)
                    z_base = radius * math.sin(theta) - 15.0 
                    y_base = height
                elif shape == "Nanopillar":
                    radius = 8.0
                    theta = (col / grid_dim) * 2 * math.pi
                    height = (row * spacing) - 10.0
                    x_base = radius * math.cos(theta)
                    y_base = radius * math.sin(theta)
                    z_base = height - 20.0 

                for atom in atoms:
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    final_pos = np.array([pos.x, pos.y, pos.z]) + np.array([x_base, y_base, z_base])
                    
                    line = (f"ATOM  {atom_serial:>5}  {atom.GetSymbol():<3} "
                            f"SUR B {res_serial:>3}    "
                            f"{final_pos[0]:>8.3f}{final_pos[1]:>8.3f}{final_pos[2]:>8.3f}"
                            f"  1.00  0.00           {atom.GetSymbol():>2}")
                    pdb_lines.append(line)
                    atom_serial += 1
                res_serial += 1
        return "\n".join(pdb_lines)