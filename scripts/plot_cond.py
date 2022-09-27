from concurrent.futures import process
import os, requests, argparse
from urllib.parse import quote_plus
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_URL = "http://localhost:3030/paperDB?query="
PREFIXES = """
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

def chem_tranformation(row, info):
    chem1 = list(info).index('GENESIS_SMPL_017-CHEBI_32588-UO_0000063')
    chem2 = list(info).index('GENESIS_SMPL_017-CHEBI_29372-UO_1000165')
    if not isinstance(row[chem1], str):
        return 1
    elif not isinstance(row[chem2], str):
        return 2
    else:
        return 0

def plot_conditions(args):

    dataset_query = ""
    proc_query = ""
    if len(args.dataset_labels[0]) > 0:
        data_labels = "'" + "' '".join(args.dataset_labels) + "'"
        dataset_query = f"""?dataset a obo:IAO_0000100.
            ?dataset rdfs:label ?data_labels.
            VALUES ?data_labels {{{data_labels}}}.
            ?trans obo:OBI_0000299 ?dataset.
                        """
    if len(args.process_labels[0]) > 0:
        proc_labels = "'" + "' '".join(args.process_labels) + "'"
        proc_query = f"VALUES ?label {{{proc_labels}}}."
    elif len(args.ignore_process[0]) > 0:
        ignored_labels = "'" + "' '".join(args.ignore_process) + "'"
        proc_query = f"""MINUS {{
            VALUES ?label {{{ignored_labels}}}.
            }}."""

    sel = f"""
        SELECT DISTINCT ?label ?spec ?num_val ?unit
                            ?growth_medium ?chem_name 
                                (COUNT(?sampling) AS ?exp_count) WHERE {{
            {dataset_query}
            ?trans a tso:GENESIS_TRANSCRIPTOMETRY_001.
            ?trans obo:OBI_0000293 ?sp.
            ?sampling obo:OBI_0000299 ?sp.
            ?sampling ccpo:GENESIS_SMPL_503 ?reg.
            ?proc ccpo:GENESIS_SMPL_508 ?reg.
            ?proc rdfs:label ?label.
            {proc_query}
            ?reg ccpo:GENESIS_SMPL_512 ?growth_medium.
            ?reg ?has_what ?what.
            ?what a ?spec.
            ?what obo:OBI_0001938 ?vs.
            ?vs obo:OBI_0001937 ?num_val.
            ?vs obo:IAO_0000039 ?unit.
            OPTIONAL {{
                ?what a ccpo:GENESIS_SMPL_017.
                ?what obo:IAO_0000136 ?chem.
                ?chem a ?chem_name
            }}
        }} GROUP BY ?label ?spec ?num_val ?unit ?growth_medium 
                            ?chem_name ?exp_count
        """

    query = BASE_URL + quote_plus(PREFIXES + sel)
    r = requests.get(query)

    experiments = {}
    if r.ok:
        # extract results from json response into dicts with keys with
        # information about condition type and unit, probably only suitable
        # for this plotting example
        for res in r.json()['results']['bindings']:
            cond_label = res['label']['value']
            if cond_label not in experiments.keys():
                experiments[cond_label] = {}

            experiments[cond_label]['growth_medium'] = res['growth_medium']\
                ['value'].split('#')[-1]
            cond_type = res['spec']['value'].split('/')[-1].split('#')[-1]
            unit = res['unit']['value'].split('/')[-1].split('#')[-1]
            if 'chem_name' in res.keys():
                chem = res['chem_name']['value'].split('/')[-1]
                key = f"{cond_type}-{chem}-{unit}"
            else:
                key = f"{cond_type}-{unit}"
            
            experiments[cond_label][key] = float(res['num_val']['value'])
    else:
        print('r not ok....')
        for res in r:
            print(res)
        raise requests.exceptions.RequestException('Something gone wrong when '
                                                   'querying database.')

    # convert results to dataframe
    df = pd.DataFrame(experiments).fillna('No modification').transpose()

    # add simple chemical dimension reduced from chemical modifications
    df['added_chem'] = df.apply(chem_tranformation, axis=1, info=df.columns)

    if not args.ignore_chemicals:
        # df with number of experiments for specified conditions
        pivots = df.pivot_table(columns=['GENESIS_SMPL_008-UO_0000027',
                                        'GENESIS_SMPL_022-UO_0000186',
                                        'added_chem'],
                                aggfunc='size').unstack(level=2)

        p = pivots.to_numpy()
        idx = np.argwhere(pivots.notna().to_numpy()).tolist()
        cols = pivots.columns
        rows = pivots.index

        # get experimental condition values to plot
        vals = np.array([(rows[i[0]][1], rows[i[0]][0], cols[i[1]])
                         for i in idx])
        # size of scatter plot points based on number of
        # experiments for conditions
        sizes = [400*p[i[0],i[1]] for i in idx]

        # plot as 3d scatter plot
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(vals[:,0], vals[:,1], vals[:,2], s=sizes)
        ax.set_xlabel('pH', fontsize=16)
        ax.set_ylabel(r'Temperature / $\degree$C', fontsize=16)
        ax.set_zlabel('Chemical stress', fontsize=16)
        ax.set_zticks([0, 1, 2])
        ax.view_init(elev=20, azim=140)
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    else:
        # as above, but only temperature and pH,
        # ie ignore chemical modifications
        pivots = df.pivot_table(columns=['GENESIS_SMPL_008-UO_0000027',
                                        'GENESIS_SMPL_022-UO_0000186'],
                                aggfunc='size').unstack(level=0)

        p = pivots.to_numpy()
        idx = np.argwhere(pivots.notna().to_numpy()).tolist()
        cols = pivots.columns
        rows = pivots.index

        vals = np.array([(rows[i[0]], cols[i[1]]) for i in idx])
        sizes = [400*p[i[0],i[1]] for i in idx]

        # 2d scatter plot
        ax = plt.figure().gca()
        plt.scatter(vals[:,0], vals[:,1], s=sizes)
        plt.xlim(3,6)
        plt.xlabel('pH')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.ylabel(r'Temperature $\degree$C')
    
    if args.save_fig:
        fname = 'cond.eps'
        if len(args.process_labels[0]) > 0:
            fname = f"cond-w-{'-'.join(args.process_labels)}.eps"
        elif len(args.ignore_process[0]) > 0:
            fname = f"cond-wo-{'-'.join(args.ignore_process)}.eps"
        elif len(args.dataset_labels[0]) > 0:
            fname = f"cond-wds-{'-'.join(args.dataset_labels)}.eps"
        plt.savefig(fname, format='eps')
    plt.show()

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_labels', type=str, default='',
                        help='label of the specific dataset/sample, eg SCS1. '
                             'If multiple, separate with comma (,), '
                             'eg SCS1,SCS2. If blank all datasets will be '
                             'retrieved.')
    parser.add_argument('--process_labels', type=str, default='',
                        help='label of the specific culturing processes, eg '
                             'SCS. If multiple, separate with comma (,), '
                             'eg SCS,SCP. If blank all processes will be '
                             'retrieved.')
    parser.add_argument('--ignore_process', type=str , default='',
                        help='label of specific culturing processes to '
                             'ignore, eg SCS. If multiple, separate with '
                             'comma (,), eg SCS,SCP. Ignored if '
                             'process_labels are provided.')
    parser.add_argument('--ignore_chemicals', type=str2bool, default=False)
    parser.add_argument('--save_fig', type=str2bool, default=False,
                        help='Wether to save figure to file or not.')
    args = parser.parse_args()
    args.dataset_labels = args.dataset_labels.split(',')
    args.process_labels = args.process_labels.split(',')
    args.ignore_process = args.ignore_process.split(',')
    plot_conditions(args)