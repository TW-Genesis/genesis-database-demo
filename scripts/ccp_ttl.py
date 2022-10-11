import uuid
import re
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# content of medium  as list with [ref to chemical, concentration, ref to unit]
medium = [['obo:CHEBI_62946', 5, 'obo:UO_0000175'],
          ['obo:CHEBI_63036', 3, 'obo:UO_0000175'],
          ['obo:CHEBI_31795', 0.5, 'obo:UO_0000175'],
          ['obo:CHEBI_17234', 7.5, 'obo:UO_0000175'],
          ['obo:CHEBI_53514', 1, 'obo:UO_0000175']]

conditions = [['SCS', 30, 'obo:UO_0000027', 5.5, 'obo:UO_0000196', 0.05, ''],
              ['SCT', 36, 'obo:UO_0000027', 5.5, 'obo:UO_0000196'],
              ['SCP', 36, 'obo:UO_0000027', 3.5, 'obo:UO_0000196'],
              ['KCl', 30, 'obo:UO_0000027', 5.5, 'obo:UO_0000196',
               'obo:CHEBI_32588', 600, 'obo:UO_0000063'],
              ['Ana', 30, 'obo:UO_0000027', 5.5, 'obo:UO_0000196',
               'obo:CHEBI_29372', 0, 'obo:UO_1000165']]
# [name, temp, ph, volumetric flow rate (mL/h), volume, [extra chemical things]],
# each condition is given name ccp:{uuid3(NAMESPACE_OID, condition name)
conditions = [['SCS', 30, 5.5, 50, 500, []],
              ['SCT', 36, 5.5, 50, 500, []],
              ['SCP', 30, 3.5, 50, 500, []],
              ['KCl', 30, 5.5, 50, 500, [['obo:CHEBI_32588', 600,
                                         'obo:UO_0000063']]],
              ['Ana', 30, 5.5, 50, 500, [['obo:CHEBI_29372', 0,
                                         'obo:UO_1000165']]]]
# [[sample1 name, [volume, unit], [time, unit]], [...], ...], if information unavailable - enter [], each ccps is given name ccps:{uuid3(NAMESPACE_OID, sample name)}. sample name 'abcX' is assumed to be related to condition 'abc'
samplings = [['SCS1', [], [50, 'obo:UO_0000003']],
             ['SCS2', [], [50, 'obo:UO_0000003']],
             ['SCS3', [], [50, 'obo:UO_0000003']],
             ['SCT1', [], [50, 'obo:UO_0000003']],
             ['SCT2', [], [50, 'obo:UO_0000003']],
             ['SCT3', [], [50, 'obo:UO_0000003']],
             ['SCP1', [], [50, 'obo:UO_0000003']],
             ['SCP2', [], [50, 'obo:UO_0000003']],
             ['SCP3', [], [50, 'obo:UO_0000003']],
             ['KCl1', [], [50, 'obo:UO_0000003']],
             ['KCl2', [], [50, 'obo:UO_0000003']],
             ['KCl3', [], [50, 'obo:UO_0000003']],
             ['Ana1', [], [50, 'obo:UO_0000003']],
             ['Ana2', [], [50, 'obo:UO_0000003']],
             ['Ana3', [], [50, 'obo:UO_0000003']]]

medium_ref = 'https://www.doi.org/10.1002/yea.320080703'
strain_ref = 'sgds:CEN.PK'

with open(os.path.join(BASE, 'data/ccp.ttl'), 'w') as fo:
    fo.write("""@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix owl: <http://www.w3.org/2002/07/owl#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .
@prefix obo: <http://purl.obolibrary.org/obo/> .

@prefix ccp: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/cell_culturing_process#> .
@prefix ccps: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/cell_culturing_process/sample#> .
@prefix ccpo: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/cell_culturing_process/ontology#> .

@prefix ys: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast#> .
@prefix yso: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast/ontology#> .
@prefix ysg: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast/gene#> .
@prefix ysgl: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/yeast/gene/location#> .


@prefix mb: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/metabolite#> .
@prefix ms: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/mass_spec#> .
@prefix mso: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/mass_spec/ontology#> .

@prefix ts: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/transcriptomics#> .
@prefix tsg: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/transcriptomics/gene#> .
@prefix tso: <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/transcriptomics/ontology#> .

@prefix sgds: <https://www.yeastgenome.org/strain/> .

@base <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/> .\n\n\n""")
    # specify yeast strain, ys:... a yso:GENESIS_YS_006, here also an external
    # reference to sgd.
    fo.write('ys:CEN_PK113-7D a yso:GENESIS_YS_006 .\n')
    if strain_ref is not None:
        fo.write(f'ys:CEN_PK113-7D rdfs:seeAlso {strain_ref} .\n')
    
    # specify growth medium according to above, here also external reference as
    # in paper. Do some uuid for medium as well..
    fo.write('ccp:growthMedium-SM1 a ccpo:GENESIS_SMPL_003 .\n')
    if medium_ref is not None:
        fo.write(f"ccp:growthMedium-SM1 rdfs:seeAlso <{medium_ref}> .\n")
    # set for chemicals to add to avoid adding to many duplicates
    chems_to_add = set()
    for mb in medium:
        uid = uuid.uuid3(uuid.NAMESPACE_OID, mb[0])
        # print(mb[0])
        # print(uid)
        fo.write(f'ccp:growthMedium-SM1 ccpo:GENESIS_SMPL_511 mb:{uid};\n')
        # fo.write(f'mb:{uid} a {mb[0]} .\n')
        chems_to_add.add(f'mb:{uid} a {mb[0]} .\n')
        fo.write('\tccpo:GENESIS_SMPL_510 [\n')
        fo.write('\ta ccpo:GENESIS_SMPL_017 ;\n')
        fo.write('\tobo:OBI_0001938 [\n')
        fo.write('\t\ta ccpo:GENESIS_SMPL_019 ;\n')
        fo.write(f'\t\tobo:OBI_0001937 "{mb[1]}"^^xsd:double;\n')
        fo.write(f'\t\tobo:IAO_0000039 {mb[2]}\n')
        fo.write('\t];\n')
        fo.write(f'\tobo:IAO_0000136 mb:{uid}\n')
        fo.write('] .\n')

    for cult in conditions:
        uid = uuid.uuid3(uuid.NAMESPACE_OID, cult[0])
        # a cell culturing
        fo.write(f'ccp:{uid} a ccpo:GENESIS_SMPL_005;\n')
        fo.write('\tccpo:GENESIS_SMPL_504 ys:CEN_PK113-7D;\n')
        # culture process has regime, potentially more than one
        fo.write(f'\tccpo:GENESIS_SMPL_508 _:regime-{uid};\n')
        # give ccp label from paper
        fo.write(f"\trdfs:label '{cult[0]}' .\n")
        # an experimental regime
        fo.write(f'_:regime-{uid} a ccpo:GENESIS_SMPL_006;\n')
        # regime has growth medium
        fo.write('\tccpo:GENESIS_SMPL_512 ccp:growthMedium-SM1 .\n\n')
        # add three bioreactors per condition as there's three replicates for all conditions in this example
        fo.write(f'_:bioreactor-{uid}-1 a ccpo:GENESIS_SMPL_001 .\n')
        fo.write(f'_:bioreactor-{uid}-2 a ccpo:GENESIS_SMPL_001 .\n')
        fo.write(f'_:bioreactor-{uid}-3 a ccpo:GENESIS_SMPL_001 .\n')
        fo.write(f'ccp:{uid} ccpo:GENESIS_SMPL_513 _:bioreactor-{uid}-1;\n')
        fo.write(f'\tccpo:GENESIS_SMPL_513 _:bioreactor-{uid}-2;\n')
        fo.write(f'\tccpo:GENESIS_SMPL_513 _:bioreactor-{uid}-3 .\n\n')
        # specify regime based on condition above, first temp., then ph[MISSING!!!], then flow rate
        fo.write(f'_:regime-{uid} ccpo:GENESIS_SMPL_506 [\n')
        fo.write('\t\ta ccpo:GENESIS_SMPL_008;\n')
        fo.write('\t\tobo:OBI_0001938 [\n')
        fo.write('\t\t\ta ccpo:GENESIS_SMPL_015;\n')
        fo.write(f'\t\t\tobo:OBI_0001937 "{cult[1]}"^^xsd:double;\n')
        fo.write('\t\t\tobo:IAO_0000039 obo:UO_0000027\n')
        fo.write('\t\t]\n')
        fo.write('\t];\n\n')

        # pH
        fo.write('\tccpo:GENESIS_SMPL_514 [\n')
        fo.write('\t\ta ccpo:GENESIS_SMPL_022;\n')
        fo.write('\t\tobo:OBI_0001938 [\n')
        fo.write('\t\t\ta ccpo:GENESIS_SMPL_024;\n')
        fo.write(f'\t\t\tobo:OBI_0001937 "{cult[2]}"^^xsd:double;\n')
        fo.write('\t\t\tobo:IAO_0000039 ccpo:UO_0000186\n')
        fo.write('\t\t]\n')
        fo.write('\t];\n\n')
        
        # flow rate
        fo.write('\tccpo:GENESIS_SMPL_507 [\n')
        fo.write('\t\ta ccpo:GENESIS_SMPL_010;\n')
        fo.write('\t\tobo:OBI_0001938 [\n')
        fo.write('\t\t\ta ccpo:GENESIS_SMPL_013;\n')
        fo.write(f'\t\t\tobo:OBI_0001937 "{cult[3]}"^^xsd:double;\n')
        fo.write('\t\t\tobo:IAO_0000039 ccpo:GENESIS_SMPL_1001\n')
        fo.write('\t\t];\n')
        fo.write('\t\tobo:IAO_0000136 ccp:growthMedium-SM1\n')
        fo.write('\t];\n\n')

        # volume
        fo.write('\tccpo:GENESIS_SMPL_509 [\n')
        fo.write('\t\ta obo:GENESIS_SMPL_004;\n')
        fo.write('\t\tobo:OBI_0001938 [\n')
        fo.write('\t\t\ta obo:OBI_0002139;\n')
        fo.write(f'\t\t\tobo:OBI_0001937 "{cult[4]}"^^xsd:double;\n')
        fo.write('\t\t\tobo:IAO_0000039 obo:UO_0000098\n')
        fo.write('\t\t]\n\n')
        # if len(cult[5]) > 0:
        # fo.write('];\n')
        # close_prev = False
        # chems_to_add = []
        for chem in cult[5]:
            # if close_prev:
            fo.write('\t];\n')
            close_prev = True
            chem_uid = uuid.uuid3(uuid.NAMESPACE_OID, chem[0])
            # print('extra')
            # print(mb[0])
            # print(chem_uid)
            # print(chem)
            # print()
            # fo.write()
            # chems_to_add.append(f'mb:{chem_uid} a {chem[0]} .\n')
            chems_to_add.add(f'mb:{chem_uid} a {chem[0]} .\n')
            fo.write('\tccpo:GENESIS_SMPL_510 [\n')
            fo.write('\t\ta ccpo:GENESIS_SMPL_017 ;\n')
            fo.write('\t\tobo:OBI_0001938 [\n')
            fo.write('\t\t\ta ccpo:GENESIS_SMPL_019;\n')
            fo.write(f'\t\t\tobo:OBI_0001937 "{chem[1]}"^^xsd:double;\n')
            fo.write(f'\t\t\tobo:IAO_0000039 {chem[2]}\n')
            fo.write('\t\t];\n')
            fo.write(f'\t\tobo:IAO_0000136 mb:{chem_uid}\n\n')
        fo.write('\t] .\n\n')
    fo.writelines(list(chems_to_add))

    for sample in samplings:
        sample_id = uuid.uuid3(uuid.NAMESPACE_OID, sample[0])
        cult_id = uuid.uuid3(uuid.NAMESPACE_OID, re.search('[a-z]+', sample[0], flags=re.IGNORECASE).group())
        cult_replicate = re.search('[0-9]+', sample[0]).group()
        fo.write(f'ccps:{sample_id} a obo:OBI_0000747;\n')
        fo.write(f"\trdfs:label '{sample[0]}' .\n")
        fo.write(f'_:sampling-{sample_id} a ccpo:GENESIS_SMPL_021;\n')

        fo.write(f'\tccpo:GENESIS_SMPL_501 _:bioreactor-{uid}-{cult_replicate};\n')
        # modify ontology to allow to specify which culturing process sample belongs to as well
        # fo.write(f'\tccpo:FROM_CULTURING_PROCESS ccp:{cult_id};\n')
        fo.write(f'\tccpo:GENESIS_SMPL_503 _:regime-{cult_id}')
        if len(sample[1]) > 0 or len(sample[2]) > 0:
            fo.write(';\n')
        else:
            fo.write(' .\n')
        # has output:
        fo.write(f'\tobo:OBI_0000299 ccps:{sample_id};\n')
        if len(sample[1]) > 0:
            # add volume
            fo.write('\tccpo:GENESIS_SMPL_509 [\n')
            fo.write('\t\ta ccpo:GENESIS_SMPL_004;\n')
            fo.write(f'\t\tobo:IAO_0000136 ccps:{sample_id};\n')
            fo.write('\t\tobo:OBI_0001938 [\n')
            fo.write('\t\t\ta obo:OBI_0002139;\n')
            fo.write(f'\t\t\tobo:OBI_0001937 "{sample[1][0]}"^^xsd:double;\n')
            fo.write(f'\t\t\tobo:IAO_0000039 {sample[1][1]}\n')
            fo.write('\t\t]\n')
            if len(sample[2]) > 0:
                fo.write('\t];\n')
            else:
                fo.write('\t] .\n')
        if len(sample[2]) > 0:
            fo.write('\tccpo:GENESIS_SMPL_502 [\n')
            fo.write('\t\ta obo:IAO_0000416;\n')
            fo.write(f'\t\tobo:IAO_0000136 ccps:{sample_id};\n')
            fo.write('\t\tobo:OBI_0001938 [\n')
            fo.write('\t\t\ta ccpo:GENESIS_SMPL_011;\n')
            fo.write(f'\t\t\tobo:OBI_0001937 "{sample[2][0]}"^^xsd:integer;\n')
            fo.write(f'\t\t\tobo:IAO_0000039 {sample[2][1]}\n')
            fo.write('\t\t]\n')
            fo.write('\t] .\n\n')
