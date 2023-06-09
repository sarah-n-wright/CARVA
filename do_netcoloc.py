from configparser import NoOptionError
from netcoloc import netprop_zscore
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
from os.path import exists
sys.path.append('/cellar/users/snwright/Git/rat_genetics/')

from analysis_functions import load_pcnet
from updated_netcoloc_functions import *

#TODO get format of z outputs
#TODO check loading of seeds and matrices
#TODO check functions for significance
#TODO check joining the difference z score results
#TODO update node names of pcnet?

#outdir="/cellar/users/snwright/Data/RareCommon"
#common_seeds_file=outdir+"/common_seeds_30630.txt"
#rare_seeds_file=outdir+"/rare_seeds_30630.txt"
#trait="30630"
outdir=str(sys.argv[1])
common_seeds_file=str(sys.argv[2])
rare_seeds_file=str(sys.argv[3])
trait_rare=str(sys.argv[4])
trait_common=str(sys.argv[5])

common_seeds = pd.read_csv(common_seeds_file, header=None)[0].tolist()
rare_seeds = pd.read_csv(rare_seeds_file, header=None)[0].tolist()

pc_nodes, G_PC = load_pcnet()

indiv_heats = np.load(outdir+"/individual_heats.npy")
stats = {"trait_rare": trait_rare, "trait_common": trait_common}

if len(common_seeds) > 500:
    common_seeds = common_seeds[:500]
    
if len(common_seeds) < 5:
    print("Not enough common seeds")
    stats = {**stats, 'mean_nps': 0, 'null_mean_nps': 0, 'p_mean_nps': 1,
            'size': 0, 'null_size': 0, 'p_size': 1}
else:
    if exists(outdir+"/z_common_"+trait_common+".tsv") and 1 > 2:
        z_common=pd.read_csv(outdir+"/z_common_"+trait_common+".tsv", sep="\t", index_col=0, squeeze=True, header=None)
        z_common.index.name=None
    else:
        z_common, common_heat, _ = netprop_zscore.calculate_heat_zscores(indiv_heats, pc_nodes,
                                                    dict(G_PC.degree), 
                                                    common_seeds,
                                                    num_reps=1000, alpha=0.5,
                                                    minimum_bin_size=10,
                                                    random_seed=1)

        z_common.to_csv(outdir+"/z_common_"+trait_common+".tsv", sep="\t", header=False)

    if exists(outdir+"/z_rare_"+trait_rare+".tsv") and 1 > 2:
        z_rare = pd.read_csv(outdir+"/z_rare_"+trait_rare+".tsv", sep="\t", squeeze=True, index_col=0, header=None)
        z_rare.index.name=None
    else:
        z_rare, rare_heat, _ = netprop_zscore.calculate_heat_zscores(indiv_heats, pc_nodes,
                                                    dict(G_PC.degree), 
                                                    rare_seeds,
                                                    num_reps=1000, alpha=0.5,
                                                    minimum_bin_size=10,
                                                    random_seed=1)                                    

        z_rare.to_csv(outdir+"/z_rare_"+trait_rare+".tsv", sep="\t", header=False)

    ###OR###
    z_scores=pd.DataFrame(z_common, columns=["Common"]).join(pd.DataFrame(z_rare, columns=["Rare"]))
    #z_scores["NPS"] = z_scores.apply(lambda x: x.Common * x.Rare, axis=1)


    #z_scores.to_csv(outdir+"/z_scores_"+trait_rare+"_"+trait_common+".tsv", sep='\t')
    ## calculate the significances
    observed, permuted = calculate_mean_z_score_distribution(pd.DataFrame(z_common, columns=["z"]), pd.DataFrame(z_rare, columns=["z"]), 
                                                            num_reps=1000, zero_double_negatives=False, overlap_control="bin",
                                                            seed1=common_seeds, seed2=rare_seeds)

    stats["mean_nps"] = observed
    stats["null_mean_nps"] = np.mean(permuted)
    stats["p_mean_nps"] = get_p_from_permutation_results(observed, permuted)

    observed_sz, permuted_sz = calculate_expected_overlap(pd.DataFrame(z_common, columns=["z"]), pd.DataFrame(z_rare, columns=["z"]), 
                                                        z_score_threshold=3, z1_threshold=1,
                                                        z2_threshold=1, num_reps=1000, plot=False, 
                                                        overlap_control="bin",
                                                        seed1=common_seeds, seed2=rare_seeds)

    stats["size"] = observed_sz
    stats["null_size"] = np.mean(permuted_sz)
    stats["p_size"] = get_p_from_permutation_results(observed_sz, permuted_sz)

with open(outdir+"/pilot_netcoloc_results_"+trait_common+".txt", 'a') as f:
    for i, s in enumerate(stats):
        if i == 0:
            f.write(str(stats[s]))
        else:
            f.write("\t"+str(stats[s]))
    f.write("\n")


