from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions

class MaterialGenerator:
    """
    Inverse Design Engine.
    Uses RDKit Reaction SMARTS to perform virtual synthesis on surface scaffolds.
    """
    
    def __init__(self):
        # Library of Virtual Reactions (Reaction SMARTS)
        # Format: [Reactant] >> [Product]
        self.rxn_library = {
            "hydrophobic": [
                # Friedel-Crafts Alkylation (Add Methyl to Aromatic Ring)
                {"name": "Methylated", "smarts": "[c:1]>>[c:1]C", "desc": "Added Methyl group to aromatic scaffold"},
                # Phenylation (Add Benzene)
                {"name": "Phenylated", "smarts": "[*:1]>>[*:1]c1ccccc1", "desc": "Grafted Phenyl ring"},
                # Fluorination of Aliphatic Chain
                {"name": "Fluorinated", "smarts": "[C:1]>>[C:1]F", "desc": "Fluorinated aliphatic chain"}
            ],
            "hydrophilic": [
                # Hydroxylation of Aromatic
                {"name": "Phenolic", "smarts": "[c:1]>>[c:1]O", "desc": "Created Phenol group"},
                # Hydroxylation of Aliphatic
                {"name": "Hydroxylated", "smarts": "[C:1]>>[C:1]O", "desc": "Added Hydroxyl group"},
                # Pegylation (Short chain)
                {"name": "Pegylated", "smarts": "[*:1]>>[*:1]OCCOC", "desc": "Grafted short PEG chain"}
            ],
            "positive": [
                # Amination (Aromatic)
                {"name": "Aminated (Ar)", "smarts": "[c:1]>>[c:1]N", "desc": "Added Aniline-like amine"},
                # Amination (Aliphatic)
                {"name": "Aminated (Al)", "smarts": "[C:1]>>[C:1]CN", "desc": "Added primary amine tail"},
                # Guanidination
                {"name": "Guanidinated", "smarts": "[*:1]>>[*:1]NC(=N)N", "desc": "Added Arginine-like Guanidinium"}
            ],
            "negative": [
                # Carboxylation
                {"name": "Carboxylated", "smarts": "[*:1]>>[*:1]C(=O)O", "desc": "Grafted Carboxylic Acid"},
                # Sulfonation
                {"name": "Sulfonated", "smarts": "[c:1]>>[c:1]S(=O)(=O)O", "desc": "Added Sulfonate group"}
            ]
        }

    def generate_candidates(self, base_smiles: str, goal: str):
        """
        Applies chemical reactions to the base molecule to generate candidates.
        """
        base_mol = Chem.MolFromSmiles(base_smiles)
        if not base_mol: return []
        
        # Add explicit hydrogens to ensure we replace real positions if needed
        # (Though simple SMARTS often work better on implicit)
        
        candidates = []
        modifiers = self.rxn_library.get(goal, [])
        
        seen_smiles = set()
        seen_smiles.add(Chem.MolToSmiles(base_mol, isomericSmiles=True))
        
        for mod in modifiers:
            rxn = AllChem.ReactionFromSmarts(mod['smarts'])
            
            # Run Reaction
            try:
                products = rxn.RunReactants((base_mol,))
            except Exception:
                continue
                
            # Reaction might produce multiple isomers; pick unique valid ones
            for ps in products:
                for product in ps:
                    try:
                        Chem.SanitizeMol(product)
                        smi = Chem.MolToSmiles(product, isomericSmiles=True)
                        
                        if smi not in seen_smiles:
                            candidates.append({
                                "name": f"{mod['name']} Surface",
                                "smiles": smi,
                                "reason": mod['desc']
                            })
                            seen_smiles.add(smi)
                            
                            # Limit to 1 product per reaction type to avoid cluttering UI
                            break 
                    except:
                        continue
                if len(candidates) >= 5: break # Limit total candidates
            
        return candidates