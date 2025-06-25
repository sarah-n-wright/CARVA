from network_annotation import *
import argparse
import pandas as pd
import networkx as nx
import sys
sys.path.append('/cellar/users/snwright/Git/Network_Evaluation_Tools/neteval/')
from Timer import Timer

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate global network statistics')
    parser.add_argument('--uuid', type=str, help='Network UUID',default='f1bbca31-f139-11ee-9621-005056ae23aa')
    parser.add_argument('--net_name', type=str, help='Network name', default='DIP')
    parser.add_argument('--outdir', type=str, help='Output directory', default = '/cellar/users/snwright/Data/RareCommon/outputs/')
    
    args = parser.parse_args()
    
    t = Timer()
    t.start('Load Network')
    try:
        G_nd = NDExNetwork(args.uuid, args.net_name)
    except:
        try:
            G_nd = NDExNetwork(args.uuid, args.net_name)
        except:
            G_nd = NDExNetwork(args.uuid, args.net_name)
    t.end('Load Network')
    t.start('Clustering Coefficients')
    cc = nx.clustering(G_nd.G)
    # save results
    pd.DataFrame({'clustering_coeff':cc}).to_csv(os.path.join(args.outdir, f'{args.net_name}_clustering_coefficients.csv'))
    t.end('Clustering Coefficients')
    t.print_all_times()