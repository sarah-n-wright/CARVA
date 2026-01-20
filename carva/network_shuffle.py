"""
Generate shuffled networks for null model comparison.

This script creates multiple shuffled versions of a network by randomly swapping
edges while preserving the node degree distribution. Shuffled networks are used as
null models to assess whether observed network colocalization is specific to the 
biological network and not just an emergent property of the network topology.

Shuffling method:
    - Uses edge swapping to preserve degree distribution
    - Successive shuffles build on previous shuffles for increased randomization
    - Similarity to original network is calculated and saved

Default behavior:
    - Generates 3 independent runs of 20 successive shuffles each (60 total networks)
    - Number of runs can be customized with --nruns parameter
    - Number of shuffles per run can be customized with --nshuffles parameter

Outputs:
    - <prefix>_shuffnet{i}_{j}.shuffled: Shuffled network files
    - <prefix>_shuffnet{i}_similarities.txt: Edge overlap with original network for each shuffle

Usage:
    # Default (3 runs of 20 shuffles = 60 networks)
    python network_shuffle.py network.txt --o /path/to/output --nSwaps 1000

    # Custom number of runs and shuffles (5 runs of 50 shuffles = 250 networks)
    python network_shuffle.py network.txt --o /path/to/output --nSwaps 1000 --nruns 5 --nshuffles 50

Note:
    Uses neteval library for network shuffling. Core arguments parsed by neteval.parse_arguments().
"""

from neteval.shuffle_networks import shuffle_network, parse_arguments
import neteval.data_import_export_tools as dit
import argparse
import sys
import os
import pandas as pd
import networkx as nx

def successive_shuffles(G, nSwaps, nshuffs, outpath, prefix):
    """Generate a series of shuffled networks.

    Args:
        G (networkx.Graph): Original network to shuffle.
        nSwaps (int): Number of edge swaps per shuffle iteration.
        nshuffs (int): Number of successive shuffles to generate.
        outpath (str): Directory to save output files.
        prefix (str): Prefix for output file names.

    Returns:
        None. Writes shuffled networks and similarity metrics to files.
    """
    similarity_to_original = {}
    for i in range(nshuffs):
        if i == 0:
            G_shuff = shuffle_network(G, nSwaps)
        else:
            G_shuff = shuffle_network(G_shuff.copy(), nSwaps)
        similarity_to_original[i] = len(set(G.edges()).intersection(set(G_shuff.edges())))/len(set(G.edges()))
        write_network(G_shuff, outpath, prefix, suffix=f'_{i}.shuffled')
    pd.DataFrame({'sims':similarity_to_original}).to_csv(os.path.join(outpath, prefix+'_similarities.txt'))
        
        
def write_network(G, outpath, prefix, suffix):
    """Write a network to file.

    Args:
        G (networkx.Graph): Network to write.
        outpath (str): Directory to save the file.
        prefix (str): Prefix for the file name.
        suffix (str): Suffix for the file name.

    Returns:
        None. Writes network to file using neteval data export tools.
    """
    outfile = os.path.join(outpath, prefix + suffix)
    print(outfile)
    dit.write_networkx_to_file(G, outfilepath=outfile)



if __name__ == '__main__':
    args = parse_arguments(sys.argv[1:])

    # Add custom arguments for number of runs and shuffles (not part of neteval's parse_arguments)
    parser = argparse.ArgumentParser()
    parser.add_argument('--nruns', default=3, type=int, help='Number of independent shuffle runs (default = 3)')
    parser.add_argument('--nshuffles', default=20, type=int, help='Number of successive shuffles per run (default = 20)')
    custom_args, _ = parser.parse_known_args()

    prefix = os.path.split(args.datafile)[1]
    if args.verbose:
        print(args)
        print("Analysis of", args.datafile)  
    G = dit.load_edgelist_to_networkx(args.datafile, testmode=args.testMode)
    if args.verbose:
        print("Data Loaded")
    if len(G.edges) > 0:
        # Generate independent runs of successive shuffles
        for i in range(custom_args.nruns):
            successive_shuffles(G, args.nSwaps, nshuffs=custom_args.nshuffles, outpath=args.o, prefix=prefix+f'shuffnet{i}')
            if args.verbose:
                print("Network Shuffled")

    else:
        print("NO EDGES:", args.datafile)