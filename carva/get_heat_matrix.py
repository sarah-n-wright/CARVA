from netcoloc import netprop
import numpy as np
import sys
sys.path.append('/cellar/users/snwright/Git/rat_genetics/')
from analysis_functions import load_network
import os

outdir=str(sys.argv[1])
uuid=str(sys.argv[2])
name=str(sys.argv[3])

## check the ordering of the nodes in the matrix - currently assuming it is the same and always loaded the same from NDEx

G_PC = load_network(uuid)
w_prime = netprop.get_normalized_adjacency_matrix(G_PC, conserve_heat=True, weighted=False)
np.save(os.path.join(outdir, name+"_w_prime.npy"), w_prime)
indiv_heats = netprop.get_individual_heats_matrix(w_prime, alpha=0.5)
np.save(os.path.join(outdir,name+"_individual_heats.npy"), indiv_heats)