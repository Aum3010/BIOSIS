import requests
import numpy as np
import math
import streamlit as st
from collections import Counter
from Bio import Align

class EvolutionaryEngine:
    """
    Real-time Conservation Analysis using RCSB API and Shannon Entropy.
    """
    def __init__(self):
        self.search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.fasta_url = "https://www.rcsb.org/fasta/entry/"

    def _fetch_homolog_ids(self, pdb_id, identity_cutoff=90):
        """
        Query RCSB Search API for polymer entities with >90% sequence identity.
        """
        query = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 0.1,
                    "identity_cutoff": identity_cutoff / 100.0,
                    "target": "pdb_protein_sequence",
                    "value": pdb_id
                }
            },
            "return_type": "entry",
            "request_options": {
                "return_all_hits": True
            }
        }
        
        
        try:
            r_seq = requests.get(f"{self.fasta_url}{pdb_id}", timeout=5)
            if r_seq.status_code != 200: return []            
            raw_seq = r_seq.text.split('\n')[1].strip()            
            query['query']['parameters']['value'] = raw_seq            
            r = requests.post(self.search_url, json=query, timeout=10)
            if r.status_code == 200:
                return [x['identifier'] for x in r.json().get('result_set', [])][:20]
        except Exception as e:
            print(f"Homolog Fetch Error: {e}")
            return []
        return []

    def _fetch_sequences(self, pdb_ids):
        """
        Downloads FASTA sequences for the list of PDB IDs.
        """
        sequences = []
        for pid in pdb_ids:
            try:
                r = requests.get(f"{self.fasta_url}{pid}", timeout=1)
                if r.status_code == 200:
                    # Simple FASTA parser
                    lines = r.text.split('\n')
                    if len(lines) >= 2:
                        sequences.append(lines[1].strip())
            except:
                continue
        return sequences

    def _calculate_entropy_per_position(self, query_seq, homologs):
        """
        Aligns each homolog to the query and computes Shannon Entropy per position.
        """
        # Configure Aligner
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        
        # Initialize Frequency Matrix [Length x 20AA]
        # We store list of residues seen at each query position
        position_data = [[] for _ in query_seq]
        
        for subject in homologs:
            try:
                # Align subject to query
                alignment = aligner.align(query_seq, subject)[0]
                
                aligned_query = alignment[0] 
                aligned_subj = alignment[1]  
                q_idx = 0
                for i, char_q in enumerate(aligned_query):
                    char_s = aligned_subj[i]
                    
                    if char_q != '-':
                        # This corresponds to query position q_idx
                        if char_s != '-':
                            position_data[q_idx].append(char_s)
                        q_idx += 1
            except:
                continue
                
        # Calculate Entropy
        # H = - sum(p * log2(p))
        conservation_scores = []
        
        for residues in position_data:
            if not residues: 
                conservation_scores.append(0) # No data = Variable
                continue
                
            total = len(residues)
            counts = Counter(residues)
            entropy = 0.0
            
            for aa, count in counts.items():
                p = count / total
                entropy -= p * math.log2(p)
            
            # Normalize: Max entropy for 20 AA is log2(20) ~ 4.32
            # Conservation = 1 - (Entropy / Max)
            norm_cons = max(0, 1 - (entropy / 4.32))
            conservation_scores.append(norm_cons * 100) # 0-100 scale
            
        return conservation_scores

    def get_conservation_map(self, pdb_id, query_chain_seq):
        """
        Main Pipeline:
        1. Find Homologs
        2. Fetch Sequences
        3. Calculate Entropy
        4. Return Map {Residue_Index (0-based): Score}
        """
        ids = self._fetch_homolog_ids(pdb_id)
        if len(ids) < 3: return None # Not enough data
        homolog_seqs = self._fetch_sequences(ids) 
        scores = self._calculate_entropy_per_position(query_chain_seq, homolog_seqs)        
        return {i: score for i, score in enumerate(scores)}