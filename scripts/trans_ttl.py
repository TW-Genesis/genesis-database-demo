import pandas as pd
import os, re, uuid

NW = re.compile("\W")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREFIX = """@prefix dc:       <http://purl.org/dc/elements/1.1/> .
@prefix efo:      <http://www.ebi.ac.uk/efo/> .
@prefix obo:      <http://purl.obolibrary.org/obo/> .
@prefix oboInOwl: <http://www.geneontology.org/formats/oboInOwl#> .
@prefix owl:      <http://www.w3.org/2002/07/owl#> .
@prefix owl2:     <http://www.w3.org/2006/12/owl2#> .
@prefix protege:  <http://protege.stanford.edu/plugins/owl/protege#> .
@prefix ql:       <http://semweb.mmlab.be/ns/ql#> .
@prefix rdf:      <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs:     <http://www.w3.org/2000/01/rdf-schema#> .
@prefix rml:      <http://semweb.mmlab.be/ns/rml#> .
@prefix rr:       <http://www.w3.org/ns/r2rml#> .
@prefix terms:    <http://purl.org/dc/terms/> .
@prefix xsd:      <http://www.w3.org/2001/XMLSchema#> .

@prefix sgds: <https://www.yeastgenome.org/strain/> .

@prefix genesis:  <http://project-genesis.io/ontology#> .\n\n\n"""

## Import data
data = pd.read_csv(os.path.join(BASE, 'data/sce_RNA_RAW_counts.csv'),
                   header=0, index_col=0)

# set to keep track which genes has been instantiated
added_orfs = set()
cultivar = "XXX"
fo = None
cnt = 0
# either save everything to one big .ttl file or one per experiment
same_file = False
if same_file:
    fo = open(os.path.join(BASE, f'data/transcriptomics-full.ttl'), 'w')
    fo.write(PREFIX)

for sample in data.columns:
    if not same_file:
        # start write to new file
        if sample[:3] != cultivar:
            if fo is not None:
                fo.close()
            fo = open(os.path.join(BASE, f'data/transcriptomics{cnt}.ttl'), 'w')
            cnt += 1
            fo.write(PREFIX)

            cultivar = sample[:3]

    sample_id = uuid.uuid3(uuid.NAMESPACE_OID, sample)

    # Create "Transcriptomics" individual for each PROCESS x SAMPLE, link to material sample and data set
    fo.write(f'genesis:trans-{sample_id} a obo:NCIT_C153189;\n')
    fo.write(f'\tobo:OBI_0000293 genesis:sample-{sample_id};\n')
    fo.write(f'\tobo:OBI_0000299 genesis:dataset-{sample_id} .\n')

    # create dataset individual and link it to genome structure
    fo.write(f'genesis:dataset-{sample_id} a obo:IAO_0000100;\n')
    fo.write('\tobo:IAO_0000221 [\n')
    fo.write('\t\ta genesis:GEN_000006;\n')
    fo.write('\t\tobo:RO_0000080 genesis:CEN_PK113-7D ] .\n')

    # for every gene count add datum to dataset
    for j, (orf, val) in enumerate(data[sample].items()):
        orf = NW.sub("_", orf)
        if orf not in added_orfs:
            fo.write(f'genesis:gene-{orf} a obo:OGG_0000000002 .\n')
            added_orfs.add(orf)

        fo.write(f'genesis:dataset-{sample_id} obo:RO_0002351 [\n')
        fo.write('\ta genesis:GEN_000027;\n')
        fo.write('\tobo:OBI_0001938 [\n')
        fo.write('\t\ta genesis:genesis:GEN_000026;\n')
        fo.write(f'\t\tobo:IAO_0000136 genesis:gene-{orf};\n')
        fo.write(f'\t\tobo:OBI_0001937 "{val}"^^xsd:nonNegativeInteger ] ] .\n')
fo.close()