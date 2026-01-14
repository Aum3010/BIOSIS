import json
from openai import OpenAI

class LLMFactory:
    """
    Real Interface to OpenAI GPT-4o for Bio-Design.
    """
    
    def __init__(self, api_key=None):
        self.client = None
        if api_key:
            self.client = OpenAI(api_key=api_key)

    def generate_suggestion(self, residue, res_id, failure_mode, surface_type):
        """
        Calls OpenAI API to generate a mutation strategy.
        """
        if not self.client:
            return {
                "mutation_code": "ERR",
                "target_res": "UNK",
                "llm_reasoning": "OpenAI API Key missing. Enter it in the sidebar.",
                "confidence": "0%"
            }

        system_prompt = """
        You are an expert Structural Biologist. 
        Your task is to suggest point mutations to fix protein-surface interaction issues.
        Output ONLY valid JSON in this format:
        {
            "target_residue": "LYS", 
            "reasoning": "Brief explanation...", 
            "confidence": "High/Medium/Low"
        }
        """

        user_prompt = f"""
        Context: Residue {residue}-{res_id} is facing a {surface_type} surface.
        Problem Detected: {failure_mode}.
        
        Task: Suggest the single best point mutation to resolve this (e.g., flip charge, reduce steric bulk).
        Use standard 3-letter amino acid codes (e.g., ALA, ASP, LYS).
        """

        try:
            response = self.client.chat.completions.create(
                model="gpt-4o",
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                temperature=0.3,
                response_format={"type": "json_object"}
            )
            
            # Parse Response
            content = response.choices[0].message.content
            data = json.loads(content)
            
            return {
                "mutation_code": f"{residue}-{res_id} -> {data.get('target_residue', 'ALA')}",
                "target_res": data.get('target_residue', 'ALA'),
                "llm_reasoning": data.get('reasoning', 'Optimization required.'),
                "confidence": data.get('confidence', 'Medium')
            }

        except Exception as e:
            return {
                "mutation_code": "API Error",
                "target_res": "UNK",
                "llm_reasoning": str(e),
                "confidence": "0%"
            }