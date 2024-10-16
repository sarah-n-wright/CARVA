from configparser import NoOptionError
from netcoloc import netprop_zscore
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
import os
from os.path import exists
import argparse
sys.path.append('/cellar/users/snwright/Git/rat_genetics/')

from geneset_utils import *
from network_utils import *
from updated_netcoloc_functions import *
# test if there are unidentified dependencies on updated coloc functions

#TODO get format of z outputs
#TODO check loading of seeds and matrices
#TODO check functions for significance
#TODO check joining the difference z score results
#TODO update node names of pcnet?

    #outdir="/cellar/users/snwright/Data/RareCommon"
    #common_seeds_file=outdir+"/common_seeds_30630.txt"
    #rare_seeds_file=outdir+"/rare_seeds_30630.txt"
    #trait="30630"
    
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=str, help='Path to outputs')
    parser.add_argument('--indir', type=str, help='Path to inputs')
    parser.add_argument('--trait_rare', type=str, help='Trait1 to evaluate')
    parser.add_argument('--trait_common', type=str, help='Trait2 to evaluate')
    parser.add_argument('--uuid', type=str, help='UUID of network')
    parser.add_argument('--net_name', type=str, help='Name of network')
    parser.add_argument('--overlap_control', type=str, choices=['remove', 'bin'], default='remove')
    args = parser.parse_args()

    common_seeds = load_seed_genes(args.trait_common, 'common', args.indir)
    rare_seeds = load_seed_genes(args.trait_rare, 'rare', args.indir)
    
    # check if there are enough seed to start with, if not exit
    if (len(common_seeds) < 3) or (len(rare_seeds) < 3):
        print("Not enough common/rare seeds")
        n_rare = len(rare_seeds)
        n_common = len(common_seeds)
        with open(os.path.join(args.outdir, 'not_enough_seeds.txt'), 'a') as f:
            f.write(f"{args.trait_rare}\t{n_rare}\t{args.trait_common}\t{n_common}\n")
    else: # continue with the analysis
        # load the presaved node and degree info, or create from the network as needed
        pc_nodes = load_saved_network_nodes(args.indir, args.net_name)
        pc_degree = load_saved_network_degrees(args.indir, args.net_name)
        if (pc_nodes is None) or (pc_degree is None):
            create_saved_nodes_and_degrees(args.uuid, args.indir, args.net_name, nodes=pc_nodes is None, degrees=pc_degree is None)
            pc_nodes = load_saved_network_nodes(args.indir, args.net_name)
            pc_degree = load_saved_network_degrees(args.indir, args.net_name)
            assert pc_nodes is not None
            assert pc_degree is not None
        
        common_seeds = [int(x) for x in common_seeds if int(x) in pc_nodes]
        rare_seeds = [int(x) for x in rare_seeds if int(x) in pc_nodes]
        # check if there are enough seeds present in the network, if not exit
        if (len(common_seeds) < 3) or (len(rare_seeds) < 3):
            print("Not enough common/rare seeds in network")
            n_rare = len(rare_seeds)
            n_common = len(common_seeds)    
            with open(os.path.join(args.outdir, 'not_enough_seeds.txt'), 'a') as f:
                f.write(f"{args.trait_rare}\t{n_rare}\t{args.trait_common}\t{n_common}\n")
        else:
            # load the heat matrix for network propagation
            if exists(os.path.join(args.indir, args.net_name+ "_individual_heats.npy")):
                indiv_heats = np.load(os.path.join(args.indir, args.net_name+ "_individual_heats.npy"))
            else:
                raise FileNotFoundError("Individual heats not found. Must be calculated first.")
            
            stats = {"trait_rare": args.trait_rare, "trait_common": args.trait_common, 'network':args.net_name}

            if len(common_seeds) > 500:
                common_seeds = common_seeds[:500]
                
            else:
                if exists(os.path.join(args.outdir, args.trait_common + "_z_CV.tsv")):
                    z_common=pd.read_csv(os.path.join(args.outdir, args.trait_common + "_z_CV.tsv"), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_common.index.name=None
                else:
                    z_common, common_heat, _ = netprop_zscore.calculate_heat_zscores(indiv_heats, pc_nodes,
                                                                pc_degree, 
                                                                common_seeds,
                                                                num_reps=1000, alpha=0.5,
                                                                minimum_bin_size=10)

                    z_common.to_csv(os.path.join(args.outdir, args.trait_common + "_z_CV.tsv"), sep="\t", header=False)
                    z_common=pd.read_csv(os.path.join(args.outdir, args.trait_common + "_z_CV.tsv"), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_common.index.name=None

                if exists(os.path.join(args.outdir, args.trait_rare + "_z_RV.tsv")):
                    z_rare = pd.read_csv(os.path.join(args.outdir, args.trait_rare + "_z_RV.tsv"), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_rare.index.name=None
                else:
                    z_rare, rare_heat, _ = netprop_zscore.calculate_heat_zscores(indiv_heats, pc_nodes,
                                                                pc_degree, 
                                                                rare_seeds,
                                                                num_reps=1000, alpha=0.5,
                                                                minimum_bin_size=10)                                    

                    z_rare.to_csv(os.path.join(args.outdir, args.trait_rare + "_z_RV.tsv"), sep="\t", header=False)
                    z_rare = pd.read_csv(os.path.join(args.outdir, args.trait_rare + "_z_RV.tsv"), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_rare.index.name=None

                ###OR###
                z_scores=pd.DataFrame(z_common, columns=["Common"]).join(pd.DataFrame(z_rare, columns=["Rare"]))
                #z_scores["NPS"] = z_scores.apply(lambda x: x.Common * x.Rare, axis=1)


                #z_scores.to_csv(outdir+"/z_scores_"+trait_rare+"_"+trait_common+".tsv", sep='\t')
                ## calculate the significances

                observed, permuted = calculate_mean_z_score_distribution(pd.DataFrame(z_common).rename(columns={1:'z'}), pd.DataFrame(z_rare).rename(columns={1:'z'}), 
                                                                        num_reps=1000, zero_double_negatives=False, overlap_control=args.overlap_control,
                                                                        seed1=common_seeds, seed2=rare_seeds)
                print(observed, permuted)
                stats["mean_nps"] = observed
                stats["null_mean_nps"] = np.mean(permuted)
                stats["p_mean_nps"] = get_p_from_permutation_results(observed, permuted)

                observed_sz, permuted_sz = calculate_expected_overlap(pd.DataFrame(z_common).rename(columns={1:'z'}), pd.DataFrame(z_rare).rename(columns={1:'z'}), 
                                                                    z_score_threshold=3, z1_threshold=1,
                                                                    z2_threshold=1, num_reps=1000, plot=False, 
                                                                    overlap_control=args.overlap_control,
                                                                    seed1=common_seeds, seed2=rare_seeds)

                stats["size"] = observed_sz
                stats["null_size"] = np.mean(permuted_sz)
                stats["p_size"] = get_p_from_permutation_results(observed_sz, permuted_sz)

                with open(os.path.join(args.outdir, f"pilot_netcoloc_results_{args.trait_rare}_{args.trait_common}.txt"), 'w') as f:
                    for i, s in enumerate(stats):
                        if i == 0:
                            f.write(str(stats[s]))
                        else:
                            f.write("\t"+str(stats[s]))
                    f.write("\n")


