import pandas as pd
import os, re, uuid

NW = re.compile("\W")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
## Import data
data = pd.read_csv(os.path.join(BASE, 'hlicorn/data/sce_RNA_RAW_counts.csv'),
                   header=0, index_col=0)
data.columns = pd.MultiIndex.from_arrays(
    [data.columns.str[:3], data.columns.str[3:]], names=["ccp", "replicate"]
)

cultivar_uuids = {c: uuid.uuid4() for c in data.columns.get_level_values(0).unique()}

cultivar = "XXX"
fo = None
file_cnt = 0
# either save everything to one big .ttl file or one per experiment
same_file = False
if same_file:
    fo = open(f'transcriptomics-full.ttl', 'w')

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


for i, sample in enumerate(data.columns):
    if not same_file:
        # start write to new file
        if sample[0] != cultivar:
            if fo is not None:
                fo.close()
            fo = open(f'transcriptomics0{file_cnt}.ttl', 'w')
            file_cnt += 1

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

    cultivar = sample[0]
    
    sample_id = 436399724531309234 + i
    # Create "transcriptometry" individual for each PROCESS x SAMPLE
    # e.g. `ts:436392344531307381 a tso:GENESIS_TRANSCRIPTOMETRY_001 .`
    fo.write("ts:{} a tso:GENESIS_TRANSCRIPTOMETRY_001;\n".format(sample_id))

    # Link a cell culturing process sample to each "transcriptometry" individual
    # e.g. `ts:436392344531307381 obo:OBI_0000293 ccps:38e2d39f-dd47-4ea1-8d65-c436da9987b4 .`
    # make sure ref to sample corresponds to what's used elsewhere..
    culture_sample_id = uuid.uuid3(uuid.NAMESPACE_OID, ''.join(sample))
    # print("ts:{} obo:OBI_0000293 ccps:{} .".format(sample_id, cultivar_uuids[cultivar]))
    fo.write("\tobo:OBI_0000293 ccps:{} .\n".format(culture_sample_id))

    # Create dataset for each "transcriptometry" individual
    # e.g. `_:dataset-436392344531307381 a obo:IAO_0000100 .`
    fo.write("_:dataset-{} a obo:IAO_0000100;\n".format(sample_id))
    # add sample name as rdfs:label
    fo.write('\trdfs:label "{}" .\n'.format(''.join(sample)))

    # Link dataset to sample
    # `ts:436392344531307381 obo:OBI_0000299 _:dataset-436392344531307381 .`
    fo.write("ts:{} obo:OBI_0000299 _:dataset-{} .\n".format(sample_id, sample_id))

    for j, (orf, val) in enumerate(data[sample].items()):
        # Create datum for each gene count in the dataset
        # e.g. `_:geneCountDatum-436392344531307381-### a obo:IAO_0000027 .`
        orf = NW.sub("_", orf)

        fo.write('\n')
        fo.write("_:geneCountDatum-{}-{:0>4} a obo:IAO_0000027 .\n".format(sample_id, j))

        # Link data to dataset
        # e.g. `_:dataset-436392344531307381 obo:BFO_0000051 _:geneCountDatum-436392344531307381-7112 .`
        fo.write(
            "_:dataset-{} obo:BFO_0000051 _:geneCountDatum-{}-{:0>4} .\n".format(
                sample_id, sample_id, j
            )
        )

        # Create genes for each row in sample
        # e.g. `tsg:YML031C-A a tso:GENESIS_TRANSCRIPTOMETRY_005 .`
        fo.write("tsg:{} a tso:GENESIS_TRANSCRIPTOMETRY_005 .\n".format(orf))

        # Link genes to cell culturing process
        # e.g. `ccps:38e2d39f-dd47-4ea1-8d65-c436da9987b4 tso:GENESIS_TRANSCRIPTOMETRY_505 tsg:YGR192C .`
        fo.write(
            "ccps:{} tso:GENESIS_TRANSCRIPTOMETRY_505 tsg:{} .\n".format(
                culture_sample_id, orf
            )
        )

        ## Fill in actual data points
        # e.g.
        #     tsg:YGR192C a tso:GENESIS_TRANSCRIPTOMETRY_005 .
        #     _:geneCountDatum-436392344531307381-1
        fo.write("_:geneCountDatum-{}-{:0>4}\n".format(sample_id, j))
        #         a tso:GENESIS_TRANSCRIPTOMETRY_004;
        fo.write("    a tso:GENESIS_TRANSCRIPTOMETRY_004;\n")
        #         obo:IAO_0000136 tsg:YGR192C;
        fo.write("    obo:IAO_0000136 tsg:{};\n".format(orf))
        # GENESIS_TRANSCRIPTOMETRY_501 means has genome structure, shouldn't it be GENESIS_TRANSCRIPTOMETRY_801 or something?
        fo.write(
            "\n".join(
                [
                    "    obo:OBI_0001938 [",
                    "        a tso:GENESIS_TRANSCRIPTOMETRY_003;",
                    '        tso:GENESIS_TRANSCRIPTOMETRY_501 "{}"^^xsd:integer;'.format(val),
                    "    ] .\n"
                ]
            )
        )
        #     _:geneCountVS436392344531307381-1 tso:GENESIS_TRANSCRIPTOMETRY_801 287235 .
        fo.write("    _:geneCountVS{}-{} tso:GENESIS_TRANSCRIPTOMETRY_801 {} .\n".format(sample_id, j, val))

    fo.write("\n")

fo.close()