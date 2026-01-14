import requests
import streamlit as st
import time

class LiteratureScout:
    """
    The 'Oracle' Module.
    Fetches context-aware literature without forcing sidebar placement.
    """
    
    def __init__(self):
        self.openalex_url = "https://api.openalex.org/works"
        self.crossref_url = "https://api.crossref.org/works"

    def _clean_query_term(self, term):
        if not term: return ""
        if any(c in term for c in ['=', '#', '(', ')', '[']):
            return "functionalized surface"
        return term

    def generate_search_strategy(self, protein_name, surface_name):
        p = protein_name
        s = self._clean_query_term(surface_name)
        return [
            f"{p} immobilization {s}",
            f"{p} adsorption {s} surface",
            f"enzyme immobilization {s} stability"
        ]

    @staticmethod
    @st.cache_data(ttl=3600, show_spinner=False)
    def _cached_search(query):
        results = []
        # 1. Try OpenAlex
        try:
            params = {
                "search": query,
                "per-page": 5,
                "filter": "type:journal-article,is_oa:true",
                "sort": "relevance_score:desc"
            }
            r = requests.get("https://api.openalex.org/works", params=params, timeout=2)
            if r.status_code == 200:
                for item in r.json().get("results", []):
                    results.append({
                        "title": item.get("title", "Unknown"),
                        "year": item.get("publication_year", ""),
                        "doi": item.get("doi", ""),
                        "journal": item.get("primary_location", {}).get("source", {}).get("display_name", "Journal"),
                        "source": "OpenAlex"
                    })
        except: pass

        if len(results) >= 3: return results

        # 2. Try Crossref
        try:
            params = {"query": query, "rows": 3}
            r = requests.get("https://api.crossref.org/works", params=params, timeout=3)
            if r.status_code == 200:
                for item in r.json().get("message", {}).get("items", []):
                    title = item.get("title", ["Unknown"])[0]
                    results.append({
                        "title": title,
                        "year": item.get("published-print", {}).get("date-parts", [[None]])[0][0],
                        "doi": f"https://doi.org/{item.get('DOI', '')}",
                        "journal": item.get("container-title", ["Journal"])[0],
                        "source": "Crossref"
                    })
        except: pass
            
        return results

    def render_feed(self, protein, surface_raw_name):
        """
        Renders the feed in the CURRENT container (not sidebar).
        """
        st.markdown("#### 📚 Literature Context")
        
        queries = self.generate_search_strategy(protein, surface_raw_name)
        all_hits = []
        seen_titles = set()
        
        for q in queries:
            hits = self._cached_search(q)
            if hits:
                for h in hits:
                    t_norm = h['title'].lower()[:50]
                    if t_norm not in seen_titles:
                        all_hits.append(h)
                        seen_titles.add(t_norm)
                if len(all_hits) >= 3: break
        
        if not all_hits:
            st.warning("No direct matches found in OpenAccess repositories.")
            return

        for p in all_hits[:4]:
            with st.expander(f"{p['year']} | {p['journal']}"):
                st.markdown(f"**{p['title']}**")
                if p['doi']: st.markdown(f"🔗 [Read Source]({p['doi']})")