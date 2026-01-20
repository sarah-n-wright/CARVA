"""
Calculate shortest path lengths for all pairs of nodes in a network.

This script computes the shortest path length between every pair of nodes in a
network and saves the result as a matrix. The output can be used as the paths_file
parameter in network_annotation.py's NDExNetwork class to prevent repeated calculations.

Network input options:
    - Load from NDEx using --uuid
    - Load from file using --netfile (TSV format with 'Entrez_A' and 'Entrez_B' columns)

Output:
    - <outpref>path_lengths.csv: Matrix of shortest path lengths between all node pairs

Usage:
    python find_shortest_paths.py --uuid <network-uuid> --outdir /path/to/output --outpref network_name_

"""

import pandas as pd
import networkx as nx
import os
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
        G = nx.from_pandas_edgelist(pd.read_csv(args.netfile), source='Entrez_A', target='Entrez_B')
    
    else:
        G = load_network(args.uuid)
    
    path_lengths = nx.all_pairs_shortest_path_length(G)
    pd.DataFrame.from_dict(dict(path_lengths)).to_csv(os.path.join(args.outdir, args.outpref+'path_lengths.csv'))