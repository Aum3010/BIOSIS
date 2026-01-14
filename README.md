# BIOSIS: Biomolecule Interaction On Surfaces Insight System

## A Surface-Agnostic Rational Design Workbench
### Built for the AIMS Lab Technical Challenge by Aum Pandya

![BIOSIS Banner](https://img.shields.io/badge/Status-Prototype_Ready-success) ![Python](https://img.shields.io/badge/Python-3.10-blue) ![Streamlit](https://img.shields.io/badge/Streamlit-1.32-red)

# Deployed App
👉 **Click here to launch:** [BIOSIS Live App](https://biosis.streamlit.app/)

BIOSIS is a comprehensive platform that bridges **Molecular Biophysics** (Structure), **Evolutionary Biology** (Conservation), and **Surface Science** (Cheminformatics). It allows researchers to:

1. **Visualize** how specific enzymes (like *Carbonic Anhydrase*) orient on any chemical surface.
2. **Optimize** binding orientation to protect active sites using Genetic Algorithms.
3. **Translate** computational models directly into physical robot protocols for the wet lab.

---

## 🛠️ Architecture

BIOSIS uses a modular architecture where each Python script handles a specific domain of the simulation:

| Module | Responsibility | Tech Stack |
| :--- | :--- | :--- |
| **`bio_processor.py`** | Parses PDB structure, identifies Zinc active sites, and calculates Solvent Accessible Surface Area (SASA). | `Biopython`, `NumPy` |
| **`surface_engine.py`** | Converts raw SMILES strings into 3D surface lattices and calculates physicochemical properties (LogP, TPSA). | `RDKit`, `Chem` |
| **`evolutionary_engine.py`** | Fetches homologs and calculates Shannon Entropy to flag conserved residues that must not be mutated. | `RCSB API`, `Bio.Align` |
| **`optimizer.py`** | **Genetic Algorithm (GA)** that evolves the protein's rotation to maximize binding while minimizing active-site occlusion. | `SciPy`, `NumPy` |
| **`interaction_lab.py`** | The Physics Engine that scores interactions and performs **Bayesian Weight Refinement** based on user data. | `Pandas`, `SciPy Optimize` |
| **`generator.py`** | **Inverse Design** module that simulates chemical reactions (e.g., Fluorination) to suggest better surface materials. | `RDKit Reactions` |
| **`llm_factory.py`** | Interpretive layer that uses GPT-4o to suggest biological point mutations for resolving steric/electrostatic clashes. | `OpenAI API` |
| **`literature_scout.py`** | Dynamically builds queries and fetches relevant "Immobilization" papers from Open Access repositories. | `OpenAlex`, `Crossref` |
| **`app.py`** | Main streamlit setup code. Also generates **Opentrons Python protocols** for physical experiments and renders the interactive dashboard. | `Streamlit`, `Altair` |

---

## 📦 How to Run

1. **Install Dependencies**
   ```bash
   pip install -r requirements.txt
    ```
2. **Launch the Dashboard**
   ```bash
   streamlit run app.py
    ```
3. **OpenAI Integration** To use the AI-driven mutation suggestions (llm_factory.py), enter your OpenAI API key in a .env file. 
   
## Integration of Provided Scripts
----------------------------------

1.  **Biopython Structure Parsing (Original\_biopython\_mmcif...)**
    
    *   **Usage:** Refactored into **bio\_processor.py**.
        
    *   **Modification:** Converted static file parsing into dynamic stream handling (allowing PDB ID fetching). Added **SASA (Solvent Accessible Surface Area)** calculations using ShrakeRupley to differentiate between buried structural residues and surface-active residues.
        
2.  **RCSB API & Search (Original\_rcsb\_api...)**
    
    *   **Usage:** Refactored into **evolutionary\_engine.py**.
        
    *   **Modification:** Extended the basic search functionality to perform **Homology Modeling**. Instead of just listing sequences, the new engine aligns homologs and calculates **Shannon Entropy** to flag evolutionarily conserved residues (high-risk mutation targets).
        
3.  **Sequence Alignment Logic (Original\_RCSB\_Extract...)**
    
    *   **Usage:** Integrated into the **Risk Scoring Logic**.
        
    *   **Modification:** The sequence data is now cross-referenced with structural coordinates to create a comprehensive "Mutation Safety Score" visualization in the dashboard.
        

## Key Assumptions & Limitations
--------------------------------

*   **Rigid Body Assumption:** The protein is treated largely as a rigid body during rotation optimization. While we account for side-chain clashes via "Soft-Sphere" potentials, we do not simulate full backbone flexibility (Molecular Dynamics).
    
*   **Surface Representation:** Surfaces are generated as idealized 3D lattices based on SMILES strings. Surface roughness or defects are not currently modeled.
    
*   **Implicit Solvation:** Water interactions are approximated using Hydrophobicity scales (LogP) rather than explicit water molecule simulation.
    
*   **Heuristic Optimization:** The Genetic Algorithm finds a _near-optimal_ orientation, but it does not guarantee finding the global energy minimum (a common trade-off for real-time performance).
    

## Future Improvements (With More Time)
---------------------------------------

1.  **Molecular Dynamics (MD) Integration:**
    
    *   Implement **OpenMM** to allow the protein backbone to relax and fold onto the surface, capturing "induced fit" phenomena.
        
2.  **Machine Learning Feedback Loop:**
    
    *   Train a small Random Forest regressor on the user's calibration data to predict binding efficiency for unknown surfaces, moving beyond simple physics weights.
        
3.  **Advanced Surface Topology:**
    
    *   Allow users to upload .cif files for crystalline surfaces (e.g., MOFs, Zeolites) instead of just generating organic lattices from SMILES.
        
4.  **Batch Processing:**
    
    *   Enable the upload of a CSV of 100+ candidates to run a high-throughput virtual screening.

4.  **Add Zotero Intehration:**
    
    *   To add Zotero integration to the literature scout so one can directly connect zotero and add the papers to their zotero group.
  
## Implementation Details

### Phase 1

Everything begins with the protein. Before we can simulate any interaction, we must understand the structural constraints and evolutionary importance of the biomolecule to ensure we preserve its function.

**1\. Structural Logic and Active Site Preservation**

*   **Why I created this:** The most common failure mode in sensor design is immobilizing the enzyme in a way that blocks its active site. If the catalytic center faces the material, the sensor is dead on arrival.
    
*   **How it benefits the researcher:** It automatically acts as a safety guard. The researcher does not need to manually inspect coordinates; the system guarantees that the optimization process effectively ignores any orientation that compromises the Zinc atom.
    
*   **Implementation Location:** bio\_processor.py
    
*   **Code Origin:** Refactored and modularized from Original\_biopython\_mmcif\_1v9e\_TUTORIAL.ipynb.
    
*   **How it is Implemented:** The system uses Bio.PDB.MMCIFParser to parse the structure. It iterates through atoms to find those labeled ZN. It performs a Neighbor Search using a KD-Tree to identify all residues within 5.0 Angstroms of the Zinc atom.
    
*   **Action/Output:** The tool outputs a dataframe flagging Zinc Critical residues and prevents the optimization logic from using these residues as anchors.
    

**2\. Evolutionary Conservation Analysis**

*   **Why I created this:** Not all amino acids are equal. Some are structural load-bearing walls; others are decoration. Modifying a conserved residue often destroys protein stability.
    
*   **How it benefits the researcher:** It provides a risk map for protein engineering. It tells the researcher exactly which residues can be safely mutated to improve binding and which must be left alone to prevent the protein from falling apart.
    
*   **Implementation Location:** evolutionary\_engine.py
    
*   **Code Origin:** Conceptual logic from Original\_rcsb\_api\_1v9e\_TUTORIAL.ipynb, but the mathematical implementation is new.
    
*   **How it is Implemented:** The system queries the RCSB Search API to find proteins with greater than 90 percent sequence identity to the target. Instead of simply listing sequences, the code aligns them and calculates the Shannon Entropy in bits for each position using the formula H equals negative sum of p log2 p.
    
*   **Action/Output:** The tool generates a Conservation Score from 0 to 100 for every residue. In the user interface, high conservation residues are flagged as High Risk for mutation.
    

### Phase 2

Once the protein is defined, we must define the surface. Instead of relying on a static database of materials, I built a system that can understand chemistry dynamically.

**3\. Surface Agnostic Material Processing**

*   **Why I created this:** Real research involves testing novel polymers and nanomaterials that do not exist in standard databases. A fixed dropdown menu would render the tool obsolete for advanced materials science.
    
*   **How it benefits the researcher:** It allows for infinite flexibility. The researcher can test any material they can describe chemically, from simple plastics to complex graphene derivatives, without waiting for a database update.
    
*   **Implementation Location:** surface\_engine.py using rdkit.
    
*   **Motivation:** The project requirement was to be surface agnostic, moving beyond fixed lists of materials to support any chemical structure.
    
*   **How it is Implemented:** The system accepts a raw SMILES string, such as C1=CC=CC=C1 for Benzene. It computes physicochemical descriptors using rdkit.Chem.Descriptors, including LogP for hydrophobicity, TPSA for polarity, and Hydrogen Bond Donors or Acceptors. It generates a 3D topology of the surface by creating a lattice grid of the molecule using rdkit.Chem.AllChem.EmbedMolecule.
    
*   **Action/Output:** Users can input any chemical formula. The system renders it in 3D and calculates how compatible it is with the protein surface properties.
    

### Phase 3

With both the protein and surface defined, the system enters the simulation phase. This is where we visualize the docking and use algorithms to find the best fit.

**4\. Hybrid 3D Visualization System**

*   **Why I created this:** Mathematical scores are abstract. Researchers need to verify physical reality visually to check for obvious errors like steric clashes or incorrect scaling that numbers might miss.
    
*   **How it benefits the researcher:** It provides immediate visual verification. The researcher can rotate the model to see exactly how the protein sits on the lattice and confirm that the red-colored Zinc neighbors are facing upward, away from the surface.
    
*   **Implementation Location:** visualizer.py and app.py.
    
*   **Motivation:** To visually inspect steric clashes and the spatial relationship between the Zinc pocket and the material.
    
*   **How it is Implemented:** The system uses py3Dmol and stmol to render two distinct PDB blocks in the same scene: the Protein fetched from RCSB and the Surface Grid generated procedurally from the SMILES string.
    
*   **Action/Output:** An interactive 3D canvas allows the user to zoom, rotate, and visually verify if the Zinc Neighbors are facing away from the surface grid.
    

**5\. Genetic Algorithm for Orientation Optimization**

*   **Why I created this:** A protein has infinite possible orientations in 3D space. Manually rotating it to find the best fit is impossible. Random sampling is inefficient and often results in the active site being blocked.
    
*   **How it benefits the researcher:** It automates the docking process. The algorithm explores millions of orientations to find the mathematical maximum for stability while rigorously penalizing any orientation that blocks the enzyme functionality.
    
*   **Implementation Location:** optimizer.py.
    
*   **Motivation:** Knowing the protein structure is insufficient; the system must determine how it orients on the material. Random placement might block the active site.
    
*   **How it is Implemented:** A Genetic Algorithm evolves the rotation matrix using Euler angles. The fitness function rewards maximizing the contact area of Anchor Residues (residues matching the surface chemistry) and applies a large negative penalty if the Zinc active site coordinates fall within the contact zone (defined as Z less than 5.0 Angstroms).
    
*   **Action/Output:** When the user initiates optimization, the protein rotates in the 3D viewer to find the optimal binding angle that preserves enzymatic activity.
    

**6\. Virtual Linker Library**

*   **Why I created this:** Direct contact between a protein and a solid surface often forces the protein to flatten, causing it to lose its shape and function.
    
*   **How it benefits the researcher:** It introduces engineering flexibility. If direct binding fails, the researcher can instantly test if adding a PEG spacer or Alkane chain saves the enzyme, simulating a more complex functionalization strategy.
    
*   **Implementation Location:** bio\_processor.py and app.py.
    
*   **Motivation:** Direct adsorption often deactivates enzymes due to lack of flexibility. Linkers allow the protein to float above the surface.
    
*   **How it is Implemented:** Users select from a dropdown containing options like PEG-4, PEG-12, or Alkane Chain. The system physically shifts the Z coordinate of the protein away from the surface grid. It adds an Entropic Penalty to the score (longer linkers reduce stability) but removes the Steric Clash penalty from the surface.
    
*   **Action/Output:** The 3D Viewer updates to show the protein suspended at the correct distance, and the Interaction Score recalculates to favor flexible linkers for bulky proteins.
    

### Phase 4

Once a fit is found, we must stress-test it against real-world physical constraints like denaturation and crowding, which simple docking scripts ignore.

**7\. Surface Induced Denaturation Assessment**

*   **Why I created this:** A protein might stick to a surface very well (high affinity) but destroy itself in the process. Highly hydrophobic surfaces pull the protein apart, turning it inside out.
    
*   **How it benefits the researcher:** It reduces false positives. It warns the researcher that although the binding score is high, the surface chemistry is too aggressive for this specific protein structure, preventing wasted experiments.
    
*   **Implementation Location:** interaction\_lab.py.
    
*   **Motivation:** A protein might bind with high affinity but denature due to surface interactions. If a surface is too hydrophobic, the internal hydrophobic core of the protein may invert to touch the surface.
    
*   **How it is Implemented:** The system calculates the Hydrophobic Moment of the protein surface. It applies a conditional penalty logic: if the Surface LogP is high and the Protein Hydrophobicity Score is high, it flags an unfolding risk. It treats the protein as a soft body rather than a rigid structure.
    
*   **Action/Output:** A Stability Warning appears in the interface preventing false positives where sticky but destructive surfaces would otherwise score high.
    

**8\. High Density Lattice Packing Simulation**

*   **Why I created this:** Industrial sensors need to pack billions of molecules onto a chip. A large, sprawling protein takes up too much space, lowering the total signal of the device.
    
*   **How it benefits the researcher:** It optimizes for efficiency, not just binding. It allows the researcher to choose a slightly weaker binder if it has a smaller footprint, ultimately yielding a higher performance sensor.
    
*   **Implementation Location:** interaction\_lab.py.
    
*   **Motivation:** Efficiency in sensors depends on how many molecules can be packed before steric hindrance occurs.
    
*   **How it is Implemented:** The code projects the 3D Van der Waals radius of the protein onto a 2D plane using Convex Hull logic on the atomic coordinates. It divides the surface area by this footprint area, accounting for a Steric Void factor based on hexagonal packing efficiency.
    
*   **Action/Output:** The system provides a Maximum Theoretical Density metric (e.g., pmol per square cm), allowing the user to optimize for packing density rather than just binding strength.
    

### Phase 5

If the simulation shows poor results, the system does not just report failure; it suggests solutions. This transforms the tool from a passive viewer into an active design assistant.

**9\. Inverse Design and Chemical Reaction Simulation**

*   **Why I created this:** When a protein does not bind, the researcher usually has to guess which chemical group to add to the surface. This is often a trial-and-error process.
    
*   **How it benefits the researcher:** It generates hypotheses automatically. The system simulates chemical reactions to propose surface modifications (like fluorination or methylation) that would mathematically improve the binding score.
    
*   **Implementation Location:** generator.py.
    
*   **Motivation:** If a protein does not bind well to a surface, the researcher typically hypothesizes surface modifications. This module automates that process.
    
*   **How it is Implemented:** The system uses rdkit.Chem.rdChemReactions to define Virtual Reactions such as Methylation, Pegylation, and Fluorination. It applies these reactions to the input SMILES string to generate new derivative molecules. It then rescores these new molecules against the protein to check for improved binding.
    
*   **Action/Output:** The tool presents a list of Suggested Modifications, such as recommending Fluorination to increase hydrophobic retention.
    

**10\. Generative AI Mutation Suggestions**

*   **Why I created this:** Sometimes the surface is fine, but the protein needs to change. Interpreting complex electrostatic clashes is difficult for non-experts.
    
*   **How it benefits the researcher:** It acts as an expert consultant. It leverages a Large Language Model to interpret the specific failure mode and suggests precise biological mutations (like changing a positive charge to a negative one) to resolve the conflict.
    
*   **Implementation Location:** llm\_factory.py.
    
*   **Motivation:** To interpret complex failure modes like electrostatic repulsion and suggest biological point mutations.
    
*   **How it is Implemented:** The code constructs a structured prompt detailing the specific residue causing a clash and the surface properties. It sends this to the OpenAI API with instructions to act as a structural biologist.
    
*   **Action/Output:** The system suggests specific mutations, such as mutating Lysine to Glutamic Acid to reduce charge repulsion, accompanied by a confidence score.
    

### Phase 6

Finally, the system ensures the model improves over time and facilitates the physical execution of the experiment.

**11\. Bayesian Parameter Refinement**

*   **Why I created this:** No simulation is perfect. Physics models often drift from reality when applied to specific enzyme classes.
    
*   **How it benefits the researcher:** It creates a learning loop. By feeding real lab results back into the system, the tool mathematically adjusts its internal weights, becoming more accurate for the researcher's specific workflow over time.
    
*   **Implementation Location:** interaction\_lab.py.
    
*   **Motivation:** Computational models require calibration against real world bench data.
    
*   **How it is Implemented:** The user inputs real experimental pairs consisting of a Surface Type and an Observed Score. The system defines a Loss Function based on the squared difference between Predicted and Observed scores. It runs scipy.optimize.minimize to reverse engineer the weighting variables for electrostatics, hydrophobicity, and sterics.
    
*   **Action/Output:** The internal physics weights are updated to match the specific experimental reality of the user, tailoring future predictions to their specific enzyme class.
    

**12\. Automated Literature Context Retrieval**

*   **Why I created this:** Researchers need to know if their proposed experiment has been done before to avoid redundancy.
    
*   **How it benefits the researcher:** It saves research time. It dynamically pulls relevant papers for the specific protein-surface combination being analyzed, providing immediate scientific context without leaving the dashboard.
    
*   **Implementation Location:** literature\_scout.py.
    
*   **Motivation:** To assist the researcher in finding prior work related to specific immobilization strategies without leaving the tool.
    
*   **How it is Implemented:** The system dynamically constructs search queries based on the current Protein ID and Surface Name. It queries the OpenAlex and Crossref APIs to fetch open access papers.
    
*   **Action/Output:** A real time feed of relevant academic papers appears in the dashboard, specifically filtered for immobilization and stability.
    

**13\. Robotic Protocol Generation**

*   **Why I created this:** The biggest bottleneck in research is translating a computer simulation into a physical experiment. Manual pipetting is slow and error-prone.
    
*   **How it benefits the researcher:** It offers zero-friction manufacturing. It converts the simulation parameters directly into code that a robot can understand, allowing the researcher to move from computer to wet lab in minutes.
    
*   **Implementation Location:** app.py (Function: generate\_opentrons\_protocol).
    
*   **Motivation:** To bridge the gap between computational simulation and wet lab validation.
    
*   **How it is Implemented:** The system takes the simulation parameters including pH, Buffer Type, and Protein Concentration and injects them into a Python template compliant with the Opentrons API.
    
*   **Action/Output:** The user downloads a Python script that can be directly uploaded to a liquid handling robot to physically mix the buffers and protein for the experiment.