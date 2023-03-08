import os, requests, argparse
from urllib.parse import quote_plus
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_URL = "http://localhost:3030/newDB?query="
PREFIXES = """
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>

    PREFIX genesis: <http://project-genesis.io/ontology#>
    """

def numeric_representation(row):
    if row['chem'] == 'CHEBI_32588':
        return 1
    elif row['chem'] == 'CHEBI_29372':
        return 2
    else:
        return 0

def plot_conditions(args):

    dataset_query = ""
    proc_query = ""
    if len(args.dataset_labels[0]) > 0:
        data_labels = "'" + "' '".join(args.dataset_labels) + "'"
        dataset_query = f"""?sample rdfs:label ?data_labels.
            VALUES ?data_labels {{{data_labels}}}.
                        """
    if len(args.process_labels[0]) > 0:
        proc_labels = "'" + "' '".join(args.process_labels) + "'"
        proc_query = f"VALUES ?label {{{proc_labels}}}."
    elif len(args.ignore_process[0]) > 0:
        ignored_labels = "'" + "' '".join(args.ignore_process) + "'"
        proc_query = f"""MINUS {{
            VALUES ?label {{{ignored_labels}}}.
            }}."""

    sel = f"""SELECT DISTINCT ?label ?gm ?temp ?ph ?chem ?chem_conc ?chem_conc_unit (COUNT(?sample) AS ?exp_count) WHERE {{
                ?trans a obo:NCIT_C153189;
                    obo:OBI_0000293 ?sample.
                ?sample a obo:OBI_0000747;
                    ^obo:OBI_0000299 ?sampling.
                {dataset_query}
                ?regime ^genesis:GEN_000018 / genesis:GEN_000009 ?sampling;
                    rdfs:label ?label;
                    genesis:GEN_000014 ?gm;
                    genesis:GEN_000012 / obo:OBI_0001937 ?temp;
                    genesis:GEN_000013 / obo:OBI_0001937 ?ph.
                {proc_query}
                OPTIONAL {{
                    ?regime genesis:GEN_000017 ?conc.
                    ?conc obo:IAO_0000136 / a ?chem;
                        obo:IAO_0000039 ?chem_conc_unit;
                        obo:OBI_0001937 ?chem_conc.
                }}
                    
    }} GROUP BY ?label ?gm ?temp ?ph ?chem ?chem_conc ?chem_conc_unit ?exp_count"""
    query = BASE_URL + quote_plus(PREFIXES + sel)
    r = requests.get(query)

    df = pd.DataFrame()
    df.index = ['temp', 'pH', 'gm', 'chem', 'no_exp']

    if r.ok:
        # extract results from json response into dicts with keys with
        # information about condition type and unit, probably only suitable
        # for this plotting example
        for res in r.json()['results']['bindings']:
            cond_label = res['label']['value']
            df[cond_label] = [float(res['temp']['value']),
                              float(res['ph']['value']),
                              res['gm']['value'].split('#')[-1],
                              res['chem']['value'].split('/')[-1] if 'chem' in res else 'na', int(res['exp_count']['value'])]

    else:
        print('r not ok....')
        for res in r:
            print(res)
        raise requests.exceptions.RequestException('Something gone wrong when '
                                                   'querying database.')

    df = df.transpose()
    df['numeric_chemical'] = df.apply(numeric_representation, axis=1)
    pivots = df.pivot_table(columns=['temp', 'pH', 'numeric_chemical'],
                            aggfunc='size').unstack(level=2)

    if not args.ignore_chemicals:
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
        pivots = df.pivot_table(columns=['temp', 'pH'],
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
        fname = 'cond-mod.eps'
        if len(args.process_labels[0]) > 0:
            fname = f"cond-w-{'-'.join(args.process_labels)}-m.eps"
        elif len(args.ignore_process[0]) > 0:
            fname = f"cond-wo-{'-'.join(args.ignore_process)}-m.eps"
        elif len(args.dataset_labels[0]) > 0:
            fname = f"cond-wds-{'-'.join(args.dataset_labels)}-m.eps"
        fname = os.path.join(BASE, f'figs/{fname}')
        print('Figure saved to:')
        print(fname)
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