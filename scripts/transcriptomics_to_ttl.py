import pandas as pd
import owlready2 as or2
import uuid
import re

NW = re.compile("\W")

## Import data
data = pd.read_csv(
    "/Users/alexander/workspace/tw-paper/hlicorn/data/sce_RNA_RAW_counts.csv",
    header=0,
    index_col=0,
)
data.columns = pd.MultiIndex.from_arrays(
    [data.columns.str[:3], data.columns.str[3:]], names=["ccp", "replicate"]
)
# print(data.head())

## Load ontology
onto = or2.get_ontology(
    "file://"
    + "/Users/alexander/workspace/genesis-prototype/ontologies/transcriptomics-ontology.xml.owl"
).load()

cultivar_uuids = {c: uuid.uuid4() for c in data.columns.get_level_values(0).unique()}

cultivar = "XXX"
for i, sample in enumerate(data.columns):
    
    if sample[0] != cultivar:
        print("XXX" + sample[0])

    cultivar = sample[0]
    
    sample_id = 436399724531309234 + i
    # Create "transcriptometry" individual for each PROCESS x SAMPLE
    # e.g. `ts:436392344531307381 a tso:GENESIS_TRANSCRIPTOMETRY_001 .`
    print("ts:{} a tso:GENESIS_TRANSCRIPTOMETRY_001 .".format(sample_id))

    # Link a cell culturing process to each "transcriptometry" individual
    # e.g. `ts:436392344531307381 obo:OBI_0000293 ccps:38e2d39f-dd47-4ea1-8d65-c436da9987b4 .`
    print("ts:{} obo:OBI_0000293 ccps:{} .".format(sample_id, cultivar_uuids[cultivar]))
    #     TO DO: create .ttl file for these processes

    # Create dataset for each "transcriptometry" individual
    # e.g. `_:dataset-436392344531307381 a obo:IAO_0000100 .`
    print("_:dataset-{} a obo:IAO_0000100 .".format(sample_id))

    # Link dataset to sample
    # `ts:436392344531307381 obo:OBI_0000299 _:dataset-436392344531307381 .`
    print("ts:{} obo:OBI_0000299 _:dataset-{} .".format(sample_id, sample_id))

    for j, (orf, val) in enumerate(data[sample].items()):
        # Create datum for each gene count in the dataset
        # e.g. `_:geneCountDatum-436392344531307381-### a obo:IAO_0000027 .`
        orf = NW.sub("_", orf)

        print("")
        print("_:geneCountDatum-{}-{:0>4} a obo:IAO_0000027 .".format(sample_id, j))

        # Link data to dataset
        # e.g. `_:dataset-436392344531307381 obo:BFO_0000051 _:geneCountDatum-436392344531307381-7112 .`
        print(
            "_:dataset-{} obo:BFO_0000051 _:geneCountDatum-{}-{:0>4} .".format(
                sample_id, sample_id, j
            )
        )

        # Create genes for each row in sample
        # e.g. `tsg:YML031C-A a tso:GENESIS_TRANSCRIPTOMETRY_005 .`
        print("tsg:{} a tso:GENESIS_TRANSCRIPTOMETRY_005 .".format(orf))

        # Link genes to cell culturing process
        # e.g. `ccps:38e2d39f-dd47-4ea1-8d65-c436da9987b4 tso:GENESIS_TRANSCRIPTOMETRY_505 tsg:YGR192C .`
        print(
            "ccps:{} tso:GENESIS_TRANSCRIPTOMETRY_505 tsg:{} .".format(
                cultivar_uuids[cultivar], orf
            )
        )

        ## Fill in actual data points
        # e.g.
        #     tsg:YGR192C a tso:GENESIS_TRANSCRIPTOMETRY_005 .
        #     _:geneCountDatum-436392344531307381-1
        print("_:geneCountDatum-{}-{:0>4}".format(sample_id, j))
        #         a tso:GENESIS_TRANSCRIPTOMETRY_004;
        print("    a tso:GENESIS_TRANSCRIPTOMETRY_004;")
        #         obo:IAO_0000136 tsg:YGR192C;
        print("    obo:IAO_0000136 tsg:{};".format(orf))
        #         obo:OBI_0001938 [
        #             a tso:GENESIS_TRANSCRIPTOMETRY_003;
        #             tso:GENESIS_TRANSCRIPTOMETRY_501 "287235"^^xsd:integer;
        #         ] .
        print(
            "\n".join(
                [
                    "    obo:OBI_0001938 [",
                    "        a tso:GENESIS_TRANSCRIPTOMETRY_003;",
                    '        tso:GENESIS_TRANSCRIPTOMETRY_501 "{}"^^xsd:integer;'.format(val),
                    "    ] ."
                ]
            )
        )
        #     _:geneCountVS436392344531307381-1 tso:GENESIS_TRANSCRIPTOMETRY_801 287235 .
        print("    _:geneCountVS{}-{} tso:GENESIS_TRANSCRIPTOMETRY_801 {} .".format(sample_id, j, val))

    print("")