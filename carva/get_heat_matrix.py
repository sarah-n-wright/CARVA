"""
Generate network propagation matrices and node metadata for NetColoc analysis.

This script precomputes the normalized adjacency matrix (w_prime) and individual
heat diffusion matrices for a given network. These matrices are used by CARVA's
network propagation algorithms (e.g., do_carva_netcoloc.py).

Network input options:
    - Load from NDEx using --uuid
    - Load from file using --netfile (TSV format with 'Entrez_A'/'Entrez_B' or 'Node_A'/'Node_B' columns)

Optional filtering:
    - Use --filter to threshold edges by a column value (format: column_threshold, e.g., 'score')

Outputs:
    - <name>_w_prime.npy: Normalized adjacency matrix for network propagation
    - <name>_individual_heats.npy: Precomputed heat diffusion matrix (alpha=0.5 dissipation constant)
    - <name>_nodes.txt: List of network nodes
    - <name>_degrees.txt: Node degree information

Usage:
    # Load from NDEx
    python get_heat_matrix.py --uuid <network-uuid> --outdir /path/to/output --name network_name

    # Load from file with filtering and custom alpha
    python get_heat_matrix.py --netfile network.tsv --outdir /path/to/output --name network_name --filter score_10 --alpha 0.2

Note:
    The output files are required inputs for do_carva_netcoloc.py and other CARVA network analysis scripts.
"""

from netcoloc import netprop
import numpy as np
import sys
from network_utils import load_network
import os
import pandas as pd
import networkx as nx
import argparse


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Generate network diffusion matrices and node metadata.")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--uuid", help="UUID of the network to load", required=False)
    parser.add_argument("--netfile", help="Filepath of network to load", required=False)
    parser.add_argument("--name", help="Base name for output files")
    parser.add_argument("--filter", help="Filter in the format <column>_<threshold> (e.g., score_10)", required=False, default=None)
    parser.add_argument("--alpha", default=0.5, required=False, help="Dissipation constant for network propagation. Default alpha=0.5 is recommended.")
    args = parser.parse_args()
    outdir = args.outdir
    uuid = args.uuid
    name = args.name
    # Load the network
    if uuid is not None:
        G_PC = load_network(uuid)
    elif args.netfile is not None:
        df = pd.read_csv(args.netfile, sep='\t')
        try:
            G_PC = nx.from_pandas_edgelist(df, source='Entrez_A', target='Entrez_B')
        except KeyError:
            G_PC = nx.from_pandas_edgelist(df, source='Node_A', target='Node_B')

    print(len(G_PC.edges()), len(G_PC.nodes()))

    # Apply edge filtering if specified (removes edges below threshold)
    if args.filter is not None:
        filter_col, filter_th = args.filter.split('_')
        filter_th = float(filter_th)
        G_PC = nx.Graph([(u,v,d) for u,v,d in G_PC.edges(data=True) if float(d[filter_col])>=filter_th])
        print(len(G_PC.edges()), len(G_PC.nodes()))

    # Step 1: Create normalized adjacency matrix for network propagation
    w_prime = netprop.get_normalized_adjacency_matrix(G_PC, conserve_heat=True, weighted=False)
    np.save(os.path.join(outdir, name+"_w_prime.npy"), w_prime)
    # Step 2: Generate individual heat diffusion matrix (alpha=0.5 is the dissipation constant)
    indiv_heats = netprop.get_individual_heats_matrix(w_prime, alpha=args.alpha)
    np.save(os.path.join(outdir,name+"_individual_heats.npy"), indiv_heats)
    # Step 3: Save node list and degree information for downstream analyses
    pc_nodes_df = pd.DataFrame({'node': G_PC.nodes})
    pc_nodes_df.to_csv(os.path.join(outdir, name+ "_nodes.txt"), sep='\t', header=False, index=False)
    degree_map = pd.DataFrame(G_PC.degree())
    degree_map.to_csv(os.path.join(outdir, name+ "_degrees.txt"), sep='\t', header=False, index=False)
