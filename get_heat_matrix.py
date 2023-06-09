from netcoloc import netprop
import numpy as np
import sys
sys.path.append('/cellar/users/snwright/Git/rat_genetics/')
from analysis_functions import load_pcnet

outdir=str(sys.argv[1])

pc_nodes, G_PC = load_pcnet()
w_prime = netprop.get_normalized_adjacency_matrix(G_PC, conserve_heat=True, weighted=False)
np.save(outdir+"/w_prime.npy", w_prime)
indiv_heats = netprop.get_individual_heats_matrix(w_prime, alpha=0.5)
np.save(outdir+"/individual_heats.npy", indiv_heats)