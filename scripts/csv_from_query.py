import requests
from urllib.parse import quote_plus
import pandas as pd
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cultivars = {}
# change base_url to match database server
base_url = "http://localhost:3030/newDB?query="
PREFIXES = """
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>

    PREFIX genesis: <http://project-genesis.io/ontology#>
    """

# SPARQL query to retrieve all gene counts together with the corresponding
# gene name and dataset label.
# Modify this if querying for specific experimental conditions

sparql = """
    SELECT ?data_label ?gene ?cnt WHERE {
        ?t a obo:NCIT_C153189;
            obo:OBI_0000293 / rdfs:label ?data_label;
            obo:OBI_0000299 ?dataset.
        ?dataset obo:RO_0002351 / obo:OBI_0001938 ?vs.
        ?vs obo:IAO_0000136 ?gene;
            obo:OBI_0001937 ?cnt.
    }
    """

query = base_url + quote_plus(PREFIXES + sparql)
print('running query...')
r = requests.get(query)
print(r)
if r.ok:
    gs_matrix = {}
    for entry in r.json()['results']['bindings']:
        data_label = entry['data_label']['value'].replace(' ', '')
        gene = entry['gene']['value'].split('gene-')[-1].replace('_', '-')
        cnt = int(entry['cnt']['value'])
        try:
            gs_matrix[data_label][gene] = cnt
        except KeyError:
            gs_matrix[data_label] = {gene: cnt}
    df = pd.DataFrame(gs_matrix)
    fp = os.path.join(BASE, 'data/recreated.csv')
    df.to_csv(fp)
    print('saved to csv file at:')
    print(fp)
else:
    print("Something's wrong with the http response:")
    print(r.content)