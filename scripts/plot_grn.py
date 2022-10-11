import os, argparse
import pyreadr
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# List target and regulator genes to be plotted,
# if lists are empty plots entire network
REGULATOR_GENES = []
TARGET_GENES = []

TARGET_GENES = ['YDR171W', 'YPL240C', 'YLR259C', 'YNL123W']

with open(os.path.join(BASE, 'data/regulators.txt'), 'r') as fi:
    REGULATOR_GENES = fi.read().splitlines()

def plot_grn(args):
    rfile = pyreadr.read_r(
        os.path.join(BASE, f'data/reg_frames{args.cond}.Rdata'))
    df_grn = pd.concat([rfile['activator_frame'], rfile['repressor_frame']])
    if len(REGULATOR_GENES) > 0:
        df_grn = df_grn.loc[df_grn['regulators'].isin(REGULATOR_GENES)]
    if len(TARGET_GENES) > 0:
        df_grn = df_grn.loc[df_grn['targets'].isin(TARGET_GENES)]
    print('activators: ', len(df_grn.loc[df_grn['type'] == 'activator'].index))
    print('reps: ', len(df_grn.loc[df_grn['type'] == 'repressor'].index))

    g = ig.Graph.TupleList(df_grn.itertuples(index=False), directed=True,
                           edge_attrs="type")

    if args.add_unconnected:
        for tg in TARGET_GENES:
            if tg not in g.vs['name']:
                g.add_vertices(1)
                g.vs[-1]['name'] = tg

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

    # get indices for targets, for sizes when plotting
    targets = set()
    for vs in g.vs:
        if len(targets) == len(TARGET_GENES):
            break
        if vs['name'] in TARGET_GENES:
            targets.add(vs.index)

    if len(g.vs) > 15:
        label_size = [9] * len(g.vs)
        vertex_size = [15] * len(g.vs)
        margins = 35
        for i in targets:
            label_size[i] = 15
            vertex_size[i] = 30
    else:
        label_size = [15] * len(g.vs)
        vertex_size = [23] * len(g.vs)
        margins = 70
        for i in targets:
            label_size[i] = 20
            vertex_size[i] = 37
    
    g.vs["label"] = g.vs["name"]
    layout = g.layout(args.layout)
    if args.save_fig:
        ig.plot(g, layout=layout, vertex_label_size=label_size,
                target=os.path.join(BASE,
                                    f'figs/grn{args.cond}-{args.layout}.pdf'), 
                vertex_color=v_cols, vertex_size=vertex_size,
                vertex_label_dist=1.3, edge_width=2.5, margin=margins,
                edge_color=[cols[t] for t in g.es['type']])
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
                             'leave blank to plot grn for all conditions. '
                             'Make sure grn for specified conditions exists.')
    parser.add_argument('--layout', default='auto', type=str,
                        help='Layoyt of plotted graph, eg kamada_kawai, see '
                             'https://igraph.org/python/tutorial/latest/'
                             'tutorial.html#layouts-and-plotting '
                             'for more instructions')
    parser.add_argument('--save_fig', default=True, type=str2bool,
                        help='Whether graph should be saved directly to file '
                             'or showed as plt figure. Note, figures does not '
                             'look exactly the same.')
    parser.add_argument('--add_unconnected', type=str2bool, default=False,
                        help='Add unconnected targets to displayed graph.')
    
    args = parser.parse_args()
    if args.cond != '' and args.cond[0] != '-':
        args.cond = '-' + args.cond
    plot_grn(args)
