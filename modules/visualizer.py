from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
import numpy as np
import math
from Bio.PDB import PDBParser
import io

# Silence RDKit warnings
rdBase.DisableLog('rdApp.*')

class SurfaceVisualizer:
    
    def _is_uff_compatible(self, mol):
        allowed_atoms = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in allowed_atoms:
                return False
        return True

    def generate_topology(self, smiles: str, shape: str = "Flat", density: str = "Single", z_offset: float = -15.0) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return ""
        mol = Chem.AddHs(mol)
        
        try:
            res = AllChem.EmbedMolecule(mol)
            if res == -1: AllChem.EmbedMolecule(mol, useRandomCoords=True)
            if self._is_uff_compatible(mol):
                try: AllChem.UFFOptimizeMolecule(mol)
                except: pass 
        except: return ""

        conf = mol.GetConformer()
        atoms = mol.GetAtoms()
        pdb_lines = []
        atom_serial = 1
        res_serial = 1
        
        grid_dim = 12 if density == "High (Monolayer)" else 6
        spacing = 5.0 
        
        for row in range(grid_dim):
            for col in range(grid_dim):
                if shape == "Flat":
                    x_base = (row * spacing) - ((grid_dim * spacing) / 2)
                    y_base = (col * spacing) - ((grid_dim * spacing) / 2)
                    z_base = z_offset
                elif shape == "Nanopore":
                    radius = 25.0
                    theta = (col / grid_dim) * 2 * math.pi
                    height = (row * spacing) - ((grid_dim * spacing) / 2)
                    x_base = radius * math.cos(theta)
                    z_base = radius * math.sin(theta) + z_offset
                    y_base = height
                elif shape == "Nanopillar":
                    radius = 8.0
                    theta = (col / grid_dim) * 2 * math.pi
                    height = (row * spacing) - 10.0
                    x_base = radius * math.cos(theta)
                    y_base = radius * math.sin(theta)
                    z_base = height + z_offset - 5.0

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

    def generate_protein_lattice(self, pdb_string, spacing):
        """
        Clones the protein into a 3x3 grid centered at (0,0,0).
        """
        parser = PDBParser(QUIET=True)
        stream = io.StringIO(pdb_string)
        structure = parser.get_structure("unit", stream)
        atoms = list(structure.get_atoms())
        if not atoms: return ""
        
        base_coords = np.array([atom.get_coord() for atom in atoms])
        pdb_lines = []
        serial_counter = 1
        
        # 3x3 Grid Shifts
        grid_shifts = [
            (-1, 1), (0, 1), (1, 1),
            (-1, 0), (0, 0), (1, 0),
            (-1, -1), (0, -1), (1, -1)
        ]
        
        for (dx, dy) in grid_shifts:
            shift_vec = np.array([dx * spacing, dy * spacing, 0.0])
            is_center = (dx == 0 and dy == 0)
            chain_id = "A" if is_center else "B" # Use B for neighbors
            new_coords = base_coords + shift_vec
            
            for j, atom in enumerate(atoms):
                x, y, z = new_coords[j]
                # PDB Format manual write
                res_seq = atom.get_parent().id[1]
                # Ensure res_seq fits in 4 chars
                res_seq_str = str(res_seq)[-4:]
                
                line = (f"ATOM  {serial_counter:>5} {atom.name:<4}{' '}"
                        f"{atom.get_parent().resname:>3} {chain_id}{res_seq_str:>4}    "
                        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  {atom.bfactor:>5.2f} {atom.occupancy:>5.2f}           {atom.element:>2}")
                pdb_lines.append(line)
                serial_counter += 1
                
        return "\n".join(pdb_lines)