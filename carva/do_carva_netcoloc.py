from configparser import NoOptionError
import numpy as np
import pandas as pd
import os
from os.path import exists
import argparse

from geneset_utils import *
from network_utils import *
from netcoloc.netcoloc_utils import Seeds, get_degree_binning, Timer
import netcoloc.netprop_zscore as netprop_zscore
from netcoloc.netprop import *
from netcoloc.network_colocalization import *

    
def create_file_suffix(quant, transform, normalization, suff):
    """ Create a suffix for output files based on the parameters.
    Args:
        quant (bool): Whether the analysis is for quantitative traits.
        transform (str): Transformation applied to the scores.
        normalization (str): Normalization method applied to the scores.
        suff (str): Additional suffix to append.
    Returns:
        str: A string suffix that combines the parameters.
    """
    if quant:
        suffix ='_q'
    else:
        suffix = ''
    if transform is not None:
        suffix += '_'+transform
    if normalization is not None:
        suffix += '_'+normalization
    if (suff is not None):
        suffix += '_'+suff
    return suffix
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=str, help='Path to outputs')
    parser.add_argument('--indir', type=str, help='Path to inputs')
    parser.add_argument('--netdir', type=str, help='Directory containing precomputed network matrices', required=False, default=None)
    parser.add_argument('--trait_rare', type=str, help='Trait1 to evaluate')
    parser.add_argument('--trait_common', type=str, help='Trait2 to evaluate')
    parser.add_argument('--uuid', type=str, help='UUID of network')
    parser.add_argument('--net_name', type=str, help='Name of network')
    parser.add_argument('--overlap_control', type=str, choices=['remove', 'bin', 'None'], default='remove', help='Method to control for overlap between seed genes')
    parser.add_argument('--min-genes', type=int, default=3, help='Minimum number of seed genes required to run analysis')
    parser.add_argument('--quant', action='store_true', help='Whether the analysis is for quantitative traits (default: False)')
    parser.add_argument('--transform', type=str, default='neglog10', help='Transformation applied to the scores')
    parser.add_argument('--normalization', type=str, required=False, help='Normalization method applied to the scores')
    parser.add_argument('--suffix', type=str, default=None, help='Additional suffix to append')
    parser.add_argument('--binsize', type=int, default=10, help='Size of the bins for degree binning')
    parser.add_argument('--zcoloc', type=float, default=3, help='Z-score threshold for colocalization')
    parser.add_argument('--z1z2', type=float, default=1, help='Z-score threshold for trait 1 and 2')
    parser.add_argument('--stat_suffix', type=str, default=None, help='Additional suffix for statistics')
    parser.add_argument('--raresuff', type=str, default='_RV', help='Suffix for the file containing rare trait genes')
    parser.add_argument('--commonsuff', type=str, default='_CV', help='Suffix for the file containing common trait genes')
    args = parser.parse_args()

    t = Timer()
    t.start('Load seeds')
    # Load the seed genes
    common_seeds = Seeds(inputdata = os.path.join(args.indir, args.trait_common+args.commonsuff+'.txt'))
    rare_seeds = Seeds(inputdata = os.path.join(args.indir, args.trait_rare+args.raresuff + '.txt'))
    suffix = create_file_suffix(args.quant, args.transform, args.normalization, args.suffix)

    t.end('Load seeds')
    # check if there are enough seed to start with, if not exit
    if (len(common_seeds.genes) < args.min_genes) or (len(rare_seeds.genes) < args.min_genes):
        print("Not enough common/rare seeds")
        n_rare = len(rare_seeds.genes)
        n_common = len(common_seeds.genes)
        with open(os.path.join(args.outdir, 'not_enough_seeds.txt'), 'a') as f:
            f.write(f"{args.trait_rare}\t{n_rare}\t{args.trait_common}\t{n_common}\tinput\n")
    else: # continue with the analysis
        t.start('Load network')
        if args.netdir is None:
            netdir=args.indir
        else:
            netdir=args.netdir
        # filter to network:
        # load the presaved node and degree info, or create from the network as needed
        pc_nodes = load_saved_network_nodes(netdir, args.net_name)
        pc_degree = load_saved_network_degrees(netdir, args.net_name)
        if (pc_nodes is None) or (pc_degree is None):
            create_saved_nodes_and_degrees(args.uuid, netdir, args.net_name, nodes=pc_nodes is None, degrees=pc_degree is None)
            pc_nodes = load_saved_network_nodes(netdir, args.net_name)
            pc_degree = load_saved_network_degrees(netdir, args.net_name)
            assert pc_nodes is not None
            assert pc_degree is not None
        t.end('Load network')
        # subset the inputs to gene present in the network. 
        common_seeds.filter_seeds_by_network(pc_nodes)
        rare_seeds.filter_seeds_by_network(pc_nodes)
        
        
        # check that there are sufficient nodes in the network.
        if (len(common_seeds.genes) < args.min_genes) or (len(rare_seeds.genes) < args.min_genes):
            print("Not enough common/rare seeds in network")
            n_rare = len(rare_seeds.genes)
            n_common = len(common_seeds.genes)    
            with open(os.path.join(args.outdir, 'not_enough_seeds.txt'), 'a') as f:
                f.write(f"{args.trait_rare}\t{n_rare}\t{args.trait_common}\t{n_common}\tnetwork\n")
        else:
            t.start('Load heat matrix')
            # load the heat matrix for network propagation
            if exists(os.path.join(netdir, args.net_name+ "_individual_heats.npy")):
                indiv_heats = np.load(os.path.join(netdir, args.net_name+ "_individual_heats.npy"))
            else:
                raise FileNotFoundError("Individual heats not found. Must be calculated first.")
            t.end('Load heat matrix')
            # initialize statistics
            stats = {"trait_rare": args.trait_rare, "trait_common": args.trait_common, 'network':args.net_name, 'transform':args.transform, 'normalization':args.normalization}
            t.start('Common variants')

            ## Analyze the Common Genes using Q-NetColoc or B-NetColoc
            if args.quant:
                if exists(os.path.join(args.outdir, args.trait_common + f'_z{args.commonsuff}{suffix}.tsv')): # check if already calculated.
                    z_common=pd.read_csv(os.path.join(args.outdir, args.trait_common + f'_z{args.commonsuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_common.index.name=None
                else:
                #perform transformation and normalization of scores
                    common_seeds.transform_scores(method=args.transform)
                    if args.normalization in ['max', 'minmax', 'zscore', 'sum', 'log']:
                        common_seeds.normalize_scores(method=args.normalization)
                    t.start('Scored heat zscores')
                    z_common, common_heat, _ = netprop_zscore.calculate_scored_heat_zscores(indiv_heats, pc_nodes, pc_degree, common_seeds.scores, 
                                                        num_reps=1000, minimum_bin_size=args.binsize, verbose=True, normalize_heat=None, random_seed=None, Timer=t)
                    t.end('Scored heat zscores')
                    z_common.to_csv(os.path.join(args.outdir, args.trait_common + f'_z{args.commonsuff}{suffix}.tsv'), sep="\t", header=False)
                    z_common=pd.read_csv(os.path.join(args.outdir, args.trait_common + f'_z{args.commonsuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_common.index.name=None
                    
            
            else: # binarize the seed genes.
                if exists(os.path.join(args.outdir, args.trait_common + f'_z{args.commonsuff}{suffix}.tsv')): # check if already calculated.
                    z_common=pd.read_csv(os.path.join(args.outdir, args.trait_common + f'_z{args.commonsuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_common.index.name=None
                else:
                    if len(common_seeds.genes) > 500:
                        common_genes = common_seeds.get_top_ranked_genes(500, ascending=True)
                    else:
                        common_genes = common_seeds.genes
                    
                    z_common, common_heat, _ = netprop_zscore.calculate_heat_zscores(indiv_heats, pc_nodes,pc_degree, 
                                                common_genes,num_reps=1000, minimum_bin_size=args.binsize, alpha=0.5 )

                    z_common.to_csv(os.path.join(args.outdir, args.trait_common +f'_z{args.commonsuff}{suffix}.tsv'), sep="\t", header=False)
                    z_common=pd.read_csv(os.path.join(args.outdir, args.trait_common + f'_z{args.commonsuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_common.index.name=None
            t.end('Common variants')
            t.start('Rare variants')
            ## Analyze the Rare Genes using Q-NetColoc or B-NetColoc
            if args.quant:
                if exists(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv')):
                    z_rare = pd.read_csv(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_rare.index.name=None
                else:
                    rare_seeds.transform_scores(method=args.transform)
                    if args.normalization in ['max', 'minmax', 'zscore', 'sum', 'log']:
                        rare_seeds.normalize_scores(method=args.normalization)
                    t.start('Scored heat zscores')
                    z_rare, rare_heat, _ = netprop_zscore.calculate_scored_heat_zscores(indiv_heats, pc_nodes, pc_degree, rare_seeds.scores, 
                                                    num_reps=1000, minimum_bin_size=args.binsize, verbose=True, normalize_heat=None, random_seed=None, Timer=t)
                    t.end('Scored heat zscores')

                    z_rare.to_csv(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv'), sep="\t", header=False)
                    z_rare = pd.read_csv(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_rare.index.name=None
                    
            else:
                if exists(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv')):
                    z_rare = pd.read_csv(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_rare.index.name=None
                else:
                    if len(rare_seeds.genes) > 500:
                        rare_genes = rare_seeds.get_top_ranked_genes(500, ascending=True)
                    else:
                        rare_genes = rare_seeds.genes
                    z_rare, rare_heat, _ = netprop_zscore.calculate_heat_zscores(indiv_heats, pc_nodes,pc_degree, rare_genes,
                                                                num_reps=1000, alpha=0.5,minimum_bin_size=args.binsize)                                    

                    z_rare.to_csv(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv'), sep="\t", header=False)
                    z_rare = pd.read_csv(os.path.join(args.outdir, args.trait_rare + f'_z{args.raresuff}{suffix}.tsv'), sep="\t", index_col=0, header=None).squeeze('columns')
                    z_rare.index.name=None
            t.end('Rare variants')
            
            
            z_scores=pd.DataFrame(z_common, columns=["Common"]).join(pd.DataFrame(z_rare, columns=["Rare"]))

            # Calculate statistics
            t.start('Mean Z-score')
            observed, permuted = calculate_mean_z_score_distribution(pd.DataFrame(z_common).rename(columns={1:'z'}), pd.DataFrame(z_rare).rename(columns={1:'z'}), 
                                                                    num_reps=1000, zero_double_negatives=False, overlap_control=args.overlap_control,
                                                                    seed1=common_seeds.genes, seed2=rare_seeds.genes, quant=args.quant)
            print(observed, permuted)
            stats["mean_nps"] = observed
            stats["null_mean_nps"] = np.mean(permuted)
            stats["p_mean_nps"] = get_p_from_permutation_results(observed, permuted)
            t.end('Mean Z-score')
            observed_sz, permuted_sz = calculate_expected_overlap(pd.DataFrame(z_common).rename(columns={1:'z'}), pd.DataFrame(z_rare).rename(columns={1:'z'}), 
                                                                z_score_threshold=args.zcoloc, z1_threshold=args.z1z2,
                                                                z2_threshold=args.z1z2, num_reps=1000, plot=False, 
                                                                overlap_control=args.overlap_control,
                                                                seed1=common_seeds.genes, seed2=rare_seeds.genes)
            t.start('Size')
            stats["size"] = observed_sz
            stats["null_size"] = np.mean(permuted_sz)
            stats["p_size"] = get_p_from_permutation_results(observed_sz, permuted_sz)
            t.end('Size')
            if args.quant:
                stats_pref = 'qnetcoloc'
            else:
                stats_pref = 'netcoloc'
            if args.stat_suffix is None:
                stat_suffix=''
            else:
                stat_suffix=args.stat_suffix
            with open(os.path.join(args.outdir, f'{stats_pref}_{args.trait_rare}_{args.trait_common}_{suffix}{stat_suffix}.txt'), 'w') as f:
                for i, s in enumerate(stats):
                    if i == 0:
                        f.write(str(stats[s]))
                    else:
                        f.write("\t"+str(stats[s]))
                f.write("\n")
            t.print_all_times()
