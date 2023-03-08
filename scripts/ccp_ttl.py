import uuid
import re
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREFIXES = """@prefix dc:       <http://purl.org/dc/elements/1.1/> .
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
@prefix terms:    <http://purl.org/dc/terms/> .
@prefix xsd:      <http://www.w3.org/2001/XMLSchema#> .

@prefix sgds: <https://www.yeastgenome.org/strain/> .

@prefix genesis:  <http://project-genesis.io/ontology#> ."""

# [name, temp, ph, volumetric flow rate (mL/h), volume, [extra chemical things]],
conditions = [['SCS', 30, 5.5, 50, 500, []],
              ['SCT', 36, 5.5, 50, 500, []],
              ['SCP', 30, 3.5, 50, 500, []],
              ['KCl', 30, 5.5, 50, 500, [['obo:CHEBI_32588', 600,
                                         'obo:UO_0000063']]],
              ['Ana', 30, 5.5, 50, 500, [['obo:CHEBI_29372', 0,
                                         'obo:UO_1000165']]]]
# [[sample1 name, [volume, unit], [time, unit]], [...], ...], if information unavailable - enter [], sample name 'abcX' is assumed to be related to condition 'abc'
samplings = [['SCS1', [], [50, 'obo:UO_0000032']],
             ['SCS2', [], [50, 'obo:UO_0000032']],
             ['SCS3', [], [50, 'obo:UO_0000032']],
             ['SCT1', [], [50, 'obo:UO_0000032']],
             ['SCT2', [], [50, 'obo:UO_0000032']],
             ['SCT3', [], [50, 'obo:UO_0000032']],
             ['SCP1', [], [50, 'obo:UO_0000032']],
             ['SCP2', [], [50, 'obo:UO_0000032']],
             ['SCP3', [], [50, 'obo:UO_0000032']],
             ['KCl1', [], [50, 'obo:UO_0000032']],
             ['KCl2', [], [50, 'obo:UO_0000032']],
             ['KCl3', [], [50, 'obo:UO_0000032']],
             ['Ana1', [], [50, 'obo:UO_0000032']],
             ['Ana2', [], [50, 'obo:UO_0000032']],
             ['Ana3', [], [50, 'obo:UO_0000032']]]

medium_ref = 'https://www.doi.org/10.1002/yea.320080703'
strain_ref = 'sgds:CEN.PK'

with open(os.path.join(BASE, 'data/ccp.ttl'), 'w') as fo:
    fo.write(PREFIXES)
    fo.write('\n\n\n')
    # specify yeast strain, ys:... a yso:GENESIS_YS_006, here also an external
    # reference to sgd.
    # yeast - obo:FOODON_03411345
    fo.write('genesis:CEN_PK113-7D a obo:FOODON_03411345 .\n')
    if strain_ref is not None:
        fo.write(f'genesis:CEN_PK113-7D rdfs:seeAlso {strain_ref} .\n')
    
    # specify growth medium according to above, here also external reference as
    # in paper. Do some uuid for medium as well..
    fo.write('genesis:growthMedium-SM1 a obo:NCIT_C85504 .\n')
    if medium_ref is not None:
        fo.write(f"genesis:growthMedium-SM1 rdfs:seeAlso <{medium_ref}> .\n")
    fo.write('\n')

    for cult in conditions:
        uid = uuid.uuid3(uuid.NAMESPACE_OID, cult[0])
        # define regime
        fo.write(f'genesis:regime-{uid} a genesis:GEN_000003 .\n')
        # define study design, what it is about, its regimes, and its sampling times
        fo.write(f'genesis:study-{uid} a obo:OBI_0500000;\n')
        fo.write('\tgenesis:GEN_000019 genesis:CEN_PK113-7D;\n')
        fo.write('\tgenesis:GEN_000010 [\n')
        fo.write('\t\ta rdf:List;\n')
        fo.write('\t\trdf:rest rdf:nil;\n')
        fo.write(f'\t\trdf:first genesis:regime-{uid} ];\n')
        fo.write('\tgenesis:GEN_000021 [\n')
        fo.write('\t\ta rdf:List;\n')
        fo.write('\t\trdf:rest rdf:nil;\n')
        fo.write('\t\trdf:first [\n')
        fo.write('\t\t\ta genesis:GEN_000007;\n')
        fo.write('\t\t\tobo:IAO_0000039 obo:UO_0000032;\n')
        fo.write('\t\t\tobo:OBI_0001937 "50"^^xsd:double ]\n')
        fo.write('\t] .\n\n')

        # specify regime conditions
        fo.write(f'genesis:regime-{uid} rdfs:label "{cult[0]}";\n')
        # growth medium
        fo.write('\tgenesis:GEN_000014 genesis:growthMedium-SM1;\n')

        # temperature
        fo.write('\tgenesis:GEN_000012 [\n')
        fo.write('\t\ta obo:OBI_0002138;\n')
        fo.write('\t\tobo:IAO_0000039 obo:UO_0000027;\n')
        fo.write(f'\t\tobo:OBI_0001937 "{cult[1]}"^^xsd:double ];\n')

        # pH
        fo.write('\tgenesis:GEN_000013 [\n')
        fo.write('\t\ta genesis:GEN_000005;\n')
        fo.write(f'\t\tobo:OBI_0001937 "{cult[2]}"^^xsd:double ];\n')

        # flow rate
        fo.write('\tgenesis:GEN_000016 [\n')
        fo.write('\t\ta obo:PDRO_0010033;\n')
        fo.write('\t\tobo:IAO_0000039 obo:NCIT_C66962;\n')
        fo.write(f'\t\tobo:OBI_0001937 "{cult[3]}"^^xsd:double ] .\n')

        # if extra added or removed chemical
        # modelled as regime - has con. val. spec. - conc. val. spec.
        # conc. val. spec. - is about - chemical
        for c in cult[5]:
            fo.write(f'genesis:regime-{uid} genesis:GEN_000017 [\n')
            fo.write('\ta genesis:GEN_000004;\n')
            fo.write('\tobo:IAO_0000136 [\n')
            fo.write(f'\t\ta {c[0]} ];\n')
            fo.write(f'\tobo:IAO_0000039 {c[2]};\n')
            fo.write(f'\tobo:OBI_0001937 "{c[1]}"^^xsd:double ] .\n')

        fo.write('\n\n')

    # sampling and cultivation processes, 3 per regime
    for sample in samplings:
        reg_id = uuid.uuid3(uuid.NAMESPACE_OID, sample[0][:3])
        sample_id = uuid.uuid3(uuid.NAMESPACE_OID, sample[0])
        fo.write(f'genesis:well-{sample_id} a obo:NCIT_C128793 .\n')
        fo.write(f'genesis:sample-{sample_id} a obo:OBI_0000747;\n')
        fo.write(f'\trdfs:label "{sample[0]}" .\n')
        fo.write(f'genesis:sp-{sample_id} a obo:OBI_0000744;\n')
        fo.write(f'\tobo:OBI_0000293 genesis:well-{sample_id};\n')
        fo.write(f'\tobo:OBI_0000299 genesis:sample-{sample_id};\n')
        fo.write('\tgenesis:GEN_000021 [\n')
        fo.write('\t\ta genesis:GEN_000029;\n')
        fo.write('\t\tobo:OBI_0001938 [\n')
        fo.write('\t\t\ta obo:GENEPIO_0002113;\n')
        fo.write(f'\t\t\tobo:IAO_0000039 {sample[2][1]};\n')
        fo.write(f'\t\t\tobo:OBI_0001937 "{sample[2][0]}"^^xsd:double ]\n')
        fo.write('\t] .\n')
        fo.write(f'genesis:ccp-{sample_id} a obo:ENVO_01001815;\n')
        fo.write(f'\tobo:OBI_0000293 genesis:well-{sample_id};\n')
        fo.write(f'\tobo:OBI_0000293 genesis:CEN_PK113-7D;\n')
        fo.write('\tobo:BFO_0000051 [\n')
        fo.write('\t\ta rdf:List;\n')
        fo.write('\t\trdf:rest rdf:nil;\n')
        fo.write('\t\trdf:first [\n')
        fo.write('\t\t\ta obo:ENVO_01001815;\n')
        fo.write(f'\t\t\tgenesis:GEN_000018 genesis:regime-{reg_id};\n')
        fo.write('\t\t\tgenesis:GEN_000014 genesis:growthMedium-SM1;\n')
        fo.write(f'\t\t\tgenesis:GEN_000009 genesis:sp-{sample_id} ]\n')
        fo.write('\t] .\n\n')

fo.close()