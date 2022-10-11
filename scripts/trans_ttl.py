import pandas as pd
import os, re, uuid

NW = re.compile("\W")
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREFIX = """@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
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

@base <http://www.semanticweb.org/rushikeshhalle/ontologies/2021/11/> .\n\n\n"""

## Import data
data = pd.read_csv(os.path.join(BASE, 'hlicorn/data/sce_RNA_RAW_counts.csv'),
                   header=0, index_col=0)
data.columns = pd.MultiIndex.from_arrays(
    [data.columns.str[:3], data.columns.str[3:]], names=["ccp", "replicate"]
)

cultivar_uuids = {c: uuid.uuid4()
                  for c in data.columns.get_level_values(0).unique()}

cultivar = "XXX"
fo = None
file_cnt = 0
# either save everything to one big .ttl file or one per experiment
same_file = False
if same_file:
    fo = open(f'transcriptomics-full.ttl', 'w')
    fo.write(PREFIX)


for i, sample in enumerate(data.columns):
    if not same_file:
        # start write to new file
        if sample[0] != cultivar:
            if fo is not None:
                fo.close()
            fo = open(f'transcriptomics0{file_cnt}.ttl', 'w')
            file_cnt += 1

            fo.write(PREFIX)

    cultivar = sample[0]
    
    sample_id = 436399724531309234 + i
    # Create "transcriptometry" individual for each PROCESS x SAMPLE
    # e.g. `ts:436392344531307381 a tso:GENESIS_TRANSCRIPTOMETRY_001 .`
    fo.write(f"ts:{sample_id} a tso:GENESIS_TRANSCRIPTOMETRY_001;\n")

    # Link a cell culturing process sample to each "transcriptometry" individual
    # e.g. `ts:436392344531307381 obo:OBI_0000293 ccps:38e2d39f-dd47-4ea1-8d65-c436da9987b4 .`
    # make sure ref to sample corresponds to what's used elsewhere..
    culture_sample_id = uuid.uuid3(uuid.NAMESPACE_OID, ''.join(sample))
    # print("ts:{} obo:OBI_0000293 ccps:{} .".format(sample_id, cultivar_uuids[cultivar]))
    fo.write(f"\tobo:OBI_0000293 ccps:{culture_sample_id} .\n")

    # Create dataset for each "transcriptometry" individual
    # e.g. `_:dataset-436392344531307381 a obo:IAO_0000100 .`
    fo.write(f"_:dataset-{sample_id} a obo:IAO_0000100;\n")
    # add sample name as rdfs:label
    fo.write(f"""\trdfs:label "{''.join(sample)}" .\n""")

    # Link dataset to sample
    # `ts:436392344531307381 obo:OBI_0000299 _:dataset-436392344531307381 .`
    fo.write(f"ts:{sample_id} obo:OBI_0000299 _:dataset-{sample_id} .\n")

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
        fo.write(f"tsg:{orf} a tso:GENESIS_TRANSCRIPTOMETRY_005 .\n")

        # Link genes to cell culturing process
        # e.g. `ccps:38e2d39f-dd47-4ea1-8d65-c436da9987b4 tso:GENESIS_TRANSCRIPTOMETRY_505 tsg:YGR192C .`
        fo.write(f"ccps:{culture_sample_id} tso:GENESIS_TRANSCRIPTOMETRY_505 "
                 f"tsg:{orf} .\n")

        ## Fill in actual data points
        # e.g.
        #     tsg:YGR192C a tso:GENESIS_TRANSCRIPTOMETRY_005 .
        #     _:geneCountDatum-436392344531307381-1
        fo.write("_:geneCountDatum-{}-{:0>4}\n".format(sample_id, j))
        #         a tso:GENESIS_TRANSCRIPTOMETRY_004;
        fo.write("\ta tso:GENESIS_TRANSCRIPTOMETRY_004;\n")
        #         obo:IAO_0000136 tsg:YGR192C;
        fo.write(f"\tobo:IAO_0000136 tsg:{orf};\n")
        # GENESIS_TRANSCRIPTOMETRY_501 means has genome structure, shouldn't it be GENESIS_TRANSCRIPTOMETRY_801 or something?
        fo.write(
            "\n".join(
                [
                    "\tobo:OBI_0001938 [",
                    "\t\ta tso:GENESIS_TRANSCRIPTOMETRY_003;",
                    "\t\ttso:GENESIS_TRANSCRIPTOMETRY_501 "
                        f'{val}"^^xsd:integer;',
                    "\t] .\n"
                ]
            )
        )
        #     _:geneCountVS436392344531307381-1 tso:GENESIS_TRANSCRIPTOMETRY_801 287235 .
        fo.write(f"\t_:geneCountVS{sample_id}-{j} "
                 f"tso:GENESIS_TRANSCRIPTOMETRY_801 {val} .\n")

    fo.write("\n")

fo.close()