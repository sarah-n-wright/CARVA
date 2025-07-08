from netcoloc import netprop
import numpy as np
import sys
from network_utils import load_network
import os
import pandas as pd
import networkx as nx

outdir=str(sys.argv[1]) # Output directory for the heat matrix and node information
uuid=str(sys.argv[2]) # NDEx UUID of the network to load
name=str(sys.argv[3]) # Name for the output files
filt=str(sys.argv[4]).split('_') # Filter specification in the format "column_threshold", e.g., "degree_10"
filter_col = filt[0] # Column to filter on, e.g., "degree"
filter_th = int(filt[1]) # Threshold for filtering, e.g., 10

# Load the network from NDEx
G_PC = load_network(uuid)
print(len(G_PC.edges()), len(G_PC.nodes()))

# Filter the network if a filter is specified
if len(filter_col) > 0:
    G_PC = nx.Graph([(u,v,d) for u,v,d in G_PC.edges(data=True) if float(d[filter_col])>=filter_th])
    print(len(G_PC.edges()), len(G_PC.nodes()))

# Calculate the heat matrix and save it
w_prime = netprop.get_normalized_adjacency_matrix(G_PC, conserve_heat=True, weighted=False)
np.save(os.path.join(outdir, name+"_w_prime.npy"), w_prime)
indiv_heats = netprop.get_individual_heats_matrix(w_prime, alpha=0.5)
np.save(os.path.join(outdir,name+"_individual_heats.npy"), indiv_heats)

# Save additional node information
pc_nodes_df = pd.DataFrame({'node': G_PC.nodes})
pc_nodes_df.to_csv(os.path.join(outdir, name+ "_nodes.txt"), sep='\t', header=False, index=False)
degree_map = pd.DataFrame(G_PC.degree())
degree_map.to_csv(os.path.join(outdir, name+ "_degrees.txt"), sep='\t', header=False, index=False)