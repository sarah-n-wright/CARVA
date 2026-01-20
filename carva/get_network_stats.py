"""
Calculate clustering coefficients for all nodes in a network.

This script loads a network from NDEx and computes the clustering coefficient
for each node. The output can be used as the clustering_file parameter in
network_annotation.py's NDExNetwork class to prevent repeated calculations.

Network input:
    - Loads from NDEx using --uuid parameter
    - Uses NDExNetwork class from network_annotation module

Output:
    - <net_name>_clustering_coefficients.csv: CSV file with clustering coefficient for each node

Usage:
    python get_network_stats.py --uuid <network-uuid> --net_name NETWORK_NAME --outdir /path/to/output

Note:
    - The script includes retry logic (up to 3 attempts) to handle transient NDEx connection issues
    - Default values (uuid and net_name) are for a small example network for testing purposes
"""

from network_annotation import *
import argparse
import pandas as pd
import networkx as nx


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate global network statistics')
    parser.add_argument('--uuid', type=str, help='Network UUID',default='f1bbca31-f139-11ee-9621-005056ae23aa')
    parser.add_argument('--net_name', type=str, help='Network name', default='DIP')
    parser.add_argument('--outdir', type=str, help='Output directory', default = '')
    
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    # Load the network with retry logic (3 attempts) to handle connection issues
    try:
        G_nd = NDExNetwork(args.uuid, args.net_name, input_is_uuid=True)
    except:
        try:
            G_nd = NDExNetwork(args.uuid, args.net_name, input_is_uuid=True)
        except:
            G_nd = NDExNetwork(args.uuid, args.net_name, input_is_uuid=True)

    # Calculate clustering coefficient for each node in the network
    cc = nx.clustering(G_nd.G)
    # Save results
    pd.DataFrame({'clustering_coeff':cc}).to_csv(os.path.join(args.outdir, f'{args.net_name}_clustering_coefficients.csv'))
