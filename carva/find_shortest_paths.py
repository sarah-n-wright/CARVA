import pandas as pd
import networkx as nx
import os
import ndex2
import sys
import argparse

from network_utils import load_network


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=str, required=True, help='Directory to save outputs')
    parser.add_argument('--uuid', type=str, help='NDEx UUID of public network', default='d73d6357-e87b-11ee-9621-005056ae23aa')
    parser.add_argument('--outpref', type=str, help='Prefix to append to output file <outpref>path_lengths.csv', default='')
    parser.add_argument('--netfile', type=str, help='File path of network edgelist', default=None)
    args = parser.parse_args()
    
    assert (args.uuid is not None) or (args.netfile is not None), 'Either UUID or network file must be specified'
    
    if args.netfile is not None:
        G = nx.from_pandas_edgelist(args.netfile)
    
    else:
        G = load_network(args.uuid)
    
    path_lengths = nx.all_pairs_shortest_path_length(G)
    pd.DataFrame.from_dict(dict(path_lengths)).to_csv(os.path.join(outdir, args.outpref+'path_lengths.csv'))