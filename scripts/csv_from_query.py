import requests
from urllib.parse import quote_plus
import pandas as pd
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cultivars = {}
# change base_url to match database server
base_url = "http://localhost:3030/paperDB?query="
prefixes = """
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    PREFIX ccp: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/cell_culturing_process#>
    PREFIX ccps: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/cell_culturing_process/sample#>
    PREFIX ccpo: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/cell_culturing_process/ontology#>
    PREFIX ys: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast#>
    PREFIX yso: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast/ontology#>
    PREFIX ysg: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast/gene#>
    PREFIX ysgl: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast/gene/location#>
    PREFIX mb: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/metabolite#>
    PREFIX mso: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/mass_spec/ontology#>
    PREFIX ts: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/transcriptomics#>
    PREFIX tsg: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/transcriptomics/gene#>
    PREFIX tso: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/transcriptomics/ontology#>

    BASE <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/>
    """

# SPARQL query to retrieve all gene counts together with the corresponding
# gene name and dataset label.
# Modify this to query for specific experimental conditions
sel = f"""
    SELECT ?data_label ?gene ?cnt WHERE {{
        ?t rdf:type tso:GENESIS_TRANSCRIPTOMETRY_001.
        ?dataset rdf:type obo:IAO_0000100.
        ?t obo:OBI_0000299 ?dataset.
        ?dataset rdfs:label ?data_label.

        ?datum rdf:type obo:IAO_0000027.
        ?dataset obo:BFO_0000051 ?datum.
        ?datum obo:IAO_0000136 ?gene.
        ?datum obo:OBI_0001938 ?vs.
        ?vs tso:GENESIS_TRANSCRIPTOMETRY_501 ?cnt.
    }}
    """

query = base_url + quote_plus(prefixes + sel)
r = requests.get(query)

if r.ok:
    gs_matrix = {}
    for entry in r.json()['results']['bindings']:
        data_label = entry['data_label']['value'].replace(' ', '')
        gene = entry['gene']['value'].split('#')[-1].replace('_', '-')
        cnt = int(entry['cnt']['value'])
        try:
            gs_matrix[data_label][gene] = cnt
        except KeyError:
            gs_matrix[data_label] = {gene: cnt}
    df = pd.DataFrame(gs_matrix)
    fp = os.path.join(BASE, 'hlicorn/data/recreated.csv')
    df.to_csv(fp)
    print('saved to csv file at:')
    print(fp)
else:
    print("Something's wrong with the http response:")
    print(r.content)