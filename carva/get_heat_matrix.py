from netcoloc import netprop
import numpy as np
import sys
from network_utils import load_network
import os
import pandas as pd
import networkx as nx

outdir=str(sys.argv[1])
uuid=str(sys.argv[2])
name=str(sys.argv[3])
filt=str(sys.argv[4]).split('_')
filter_col = filt[0]
filter_th = int(filt[1])

## check the ordering of the nodes in the matrix - currently assuming it is the same and always loaded the same from NDEx

G_PC = load_network(uuid)
print(len(G_PC.edges()), len(G_PC.nodes()))

if len(filter_col) > 0:
    G_PC = nx.Graph([(u,v,d) for u,v,d in G_PC.edges(data=True) if float(d[filter_col])>=filter_th])
    print(len(G_PC.edges()), len(G_PC.nodes()))

w_prime = netprop.get_normalized_adjacency_matrix(G_PC, conserve_heat=True, weighted=False)
np.save(os.path.join(outdir, name+"_w_prime.npy"), w_prime)
indiv_heats = netprop.get_individual_heats_matrix(w_prime, alpha=0.5)
np.save(os.path.join(outdir,name+"_individual_heats.npy"), indiv_heats)
pc_nodes_df = pd.DataFrame({'node': G_PC.nodes})
pc_nodes_df.to_csv(os.path.join(outdir, name+ "_nodes.txt"), sep='\t', header=False, index=False)
degree_map = pd.DataFrame(G_PC.degree())
degree_map.to_csv(os.path.join(outdir, name+ "_degrees.txt"), sep='\t', header=False, index=False)