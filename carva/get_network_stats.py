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
    # Load the network, allowing for extra retries
    try:
        G_nd = NDExNetwork(args.uuid, args.net_name, input_is_uuid=True)
    except:
        try:
            G_nd = NDExNetwork(args.uuid, args.net_name, input_is_uuid=True)
        except:
            G_nd = NDExNetwork(args.uuid, args.net_name, input_is_uuid=True)

    cc = nx.clustering(G_nd.G)
    # save results
    pd.DataFrame({'clustering_coeff':cc}).to_csv(os.path.join(args.outdir, f'{args.net_name}_clustering_coefficients.csv'))
