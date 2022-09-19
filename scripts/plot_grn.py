import os
import pyreadr
import pandas as pd
import igraph as ig
import argparse
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# List regulator genes to be plotted,if list is empty plots entire network
REGULATOR_GENES = ['YBR049C', 'YLR403W']

def plot_grn(args):
    rfile = pyreadr.read_r(
        os.path.join(BASE, f'hlicorn/reg_frames{args.cond}.Rdata'))
    df_grn = pd.concat([rfile['activator_frame'], rfile['repressor_frame']])
    if len(REGULATOR_GENES) > 0:
        df_grn = df_grn.loc[df_grn['regulators'].isin(REGULATOR_GENES)]

    g = ig.Graph.TupleList(df_grn.itertuples(index=False), directed=True,
                        edge_attrs="type")

    # set colors for regulator vertices
    cols = {'activator': 'green', 'repressor': 'red'}
    v_cols = [(0.5,0.5,0.5)] * len(g.vs['name'])

    # lists to keep track of number of genes each gene
    # activates/represses, for coloring
    vs_act = [0] * len(g.vs['name'])
    vs_rep = [0] * len(g.vs['name'])
    regulators = []
    for i, edge in enumerate(g.es):
        if edge.attributes()['type'] == 'activator':
            vs_act[edge.source] += 1
            regulators.append(i)
        else:
            vs_rep[edge.source] += 1
            regulators.append(i)
    # change color of all regulators (activating/repressing at least one gene)
    # to color between red and green depending on ratio activated/repressed
    # all non regulating genes keep a gray color
    for i in range(len(vs_act)):
        act = vs_act[i]
        rep = vs_rep[i]
        if act + rep > 0:
            v_cols[i] = (float(rep) / (act + rep),
                         float(act) / (act + rep), 0.0)

    g.vs["label"] = g.vs["name"]
    layout = g.layout(args.layout)
    if args.save_fig:
        ig.plot(g, target=f'myfiletest{args.cond}.pdf', layout=layout,
                vertex_label_size=8, vertex_color=v_cols,vertex_label_dist=1.1,
                vertex_size=15, edge_color=[cols[t] for t in g.es['type']])
    else:
        fig, ax = plt.subplots()
        ig.plot(g, target=ax, layout=layout, vertex_color=v_cols,
                vertex_label_dist=1.1, vertex_size=1, vertex_label_size=5,
                edge_color=[cols[t] for t in g.es['type']])
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
    parser.add_argument('--cond', default='', type=str,
                        help='Experimental conditions to plot leave out, '
                             'leave blank to plot grn for all conditions')
    parser.add_argument('--layout', default='auto', type=str,
                        help='Layoyt of plotted graph, eg kamada_kawai, see '
                             'https://igraph.org/python/tutorial/latest/'
                             'tutorial.html#layouts-and-plotting '
                             'for more instructions')
    parser.add_argument('--save_fig', default=True, type=str2bool,
                        help='Whether graph should be saved directly to file '
                             'or showed as plt figure. Note, figures does not '
                             'look exactly the same.')
    
    args = parser.parse_args()
    if args.cond != '' and args.cond[0] != '-':
        args.cond = '-' + args.cond
    plot_grn(args)
