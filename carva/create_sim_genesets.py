"""
Generate simulated gene sets for benchmarking network colocalization.

This script creates pairs of gene sets with controlled
overlap and relevance parameters.

Key parameters:
    - overlap: Number of genes shared between the two gene sets
    - relevance: Proportion of genes from original set (1.0 = all original, 0.5 = 50% noise)
    - totalgenes: Total number of genes in each generated set
    - background: Method for selecting noise genes (currently only 'degree' matching)

Output naming convention:
    <setname>_overlap<N>_relevance<R>_totalgenes<T>_repeat<I>_background<B>_<CV|RV>.txt

Usage:
    # Qualitative mode with degree matching
    python create_sim_genesets.py --setfile genesets.txt --outdir /path/to/output \
        --netnodefile network_nodelist.txt --overlap 10 --relevance 0.8 \
        --totalgenes 100 --nrepeats 5 --background degree

    # Quantitative mode preserving scores and degree matching
    python create_sim_genesets.py --setfile trait_list.txt --outdir /path/to/output \
        --netnodefile network_nodelist.txt --overlap 10 --relevance 0.8 \
        --totalgenes 100 --nrepeats 5 --background degree --quant
"""

from configparser import NoOptionError
import numpy as np
import pandas as pd
import sys
from tqdm import tqdm
import os
from os.path import exists
import argparse
import random as rn
from collections import defaultdict
from geneset_utils import *

def partition_gene_set(all_genes, total_genes, overlap):
    """Partition a gene set into overlapping or non-overlapping subsets.
    Args:
        all_genes (list): Set of all genes to sample from.
        total_genes (int): Total number of genes to include across the two partitions.
        overlap (int): Number of genes that should be shared between the two partitions.
    Returns:
        tuple: Three lists (overlap_genes, gene_set1, gene_set2) where overlap_genes are shared
               between both sets, and gene_set1/gene_set2 are unique to each set.
    """
    gene_set = rn.sample(all_genes, total_genes)
    if overlap > 0:
        overlap_genes = rn.sample(gene_set, overlap)
        non_overlap_genes = [x for x in gene_set if x not in overlap_genes]
    else:
        overlap_genes = []
        non_overlap_genes = gene_set
    #shuffle the non_overlap_genes
    rn.shuffle(non_overlap_genes)
    #split into two even sets
    split = len(non_overlap_genes)//2
    gene_set1 = non_overlap_genes[:split]
    gene_set2 = non_overlap_genes[split:]
    return overlap_genes, gene_set1, gene_set2

def add_noise_to_gene_set(gene_set, relevance, background, all_subset_genes, netnodefile=None):
    """Replace a portion of gene set genes with degree-matched random genes.
    Args:
        gene_set (list): Original gene set to add noise to.
        relevance (float): Proportion of original genes to keep (0.0-1.0). E.g., 0.8 keeps 80% of genes and replaces 20% with noise.
        background (str): Background matching method. Currently only 'degree' is implemented.
        all_subset_genes (list): All genes already in use, to avoid selecting duplicates.
        netnodefile (str, optional): Path to network node file (used to find degree file).
    Returns:
        list: Gene set with (1-relevance) fraction of original genes replaced by degree-matched random genes.
    """
    n_to_replace = int(len(gene_set)*(1 - relevance))
    replace_genes = rn.sample(gene_set, n_to_replace)
    out_gene_set = [x for x in gene_set if x not in replace_genes]
    if background == 'degree':
        degree_file = netnodefile.replace('nodelist.txt', 'degrees.txt')
        degrees = pd.read_csv(degree_file, header=None, sep='\t', index_col=0)
        degree_bins = get_degree_bins(degrees)
        degree_bins.index.name=None
        # drop the genes that are already in the gene set
        noise_bins = degree_bins[~degree_bins.index.isin(all_subset_genes)].copy()
        noise_genes = []
        for gene in replace_genes:
            add_gene = get_matched_gene(degree_bins.loc[gene]['bin'], noise_bins)
            noise_genes.append(add_gene)
            # make sure this gene cannot be selected again
            noise_bins = noise_bins.drop(index=add_gene) 
    else:
        raise NotImplementedError('Only degree background is implemented')
    return out_gene_set + noise_genes


def add_noise_to_gene_set_quant(gene_profile, relevance, background, all_subset_genes, netnodefile=None):
    """Replace a portion of gene set genes with degree-matched random genes, preserving scores.
    Args:
        gene_profile (dict): Dictionary mapping gene IDs to scores (e.g., p-values).
        relevance (float): Proportion of original genes to keep (0.0-1.0). E.g., 0.8 keeps 80% of genes and replaces 20% with noise.
        background (str): Background matching method. Currently only 'degree' is implemented.
        all_subset_genes (list): All genes already in use, to avoid selecting duplicates.
        netnodefile (str, optional): Path to network node file (used to find degree file).
    Returns:
        dict: Gene set with (1-relevance) fraction of original genes replaced, preserving score distribution.
    """
    n_to_replace = int(len(gene_profile)*(1 - relevance))
    replace_genes = rn.sample(list(gene_profile.keys()), n_to_replace)
    replace_scores = [i for i in list(gene_profile.values())]
    out_gene_set = {x:score for x, score in gene_profile.items() if x not in replace_genes}
    if background == 'degree':
        degree_file = netnodefile.replace('nodelist.txt', 'degrees.txt')
        degrees = pd.read_csv(degree_file, header=None, sep='\t', index_col=0)
        degree_bins = get_degree_bins(degrees)
        degree_bins.index.name=None
        # drop the genes that are already in the gene set
        noise_bins = degree_bins[~degree_bins.index.isin(all_subset_genes)].copy()
        noise_genes = []
        for gene in replace_genes:
            try:
                add_gene = get_matched_gene(degree_bins.loc[gene]['bin'], noise_bins)
                noise_genes.append(add_gene)
            # make sure this gene cannot be selected again
                noise_bins = noise_bins.drop(index=add_gene)
            except KeyError:
                pass # gene not in network                
        noise_scores = {x: replace_scores[i] for i, x in enumerate(noise_genes)}
    else:
        raise NotImplementedError('Only degree background is implemented')
    return {**out_gene_set, **noise_scores}


def get_matched_gene(gene_bin, bin_data):
    """Select a random gene from the same degree bin.
    Args:
        gene_bin (int): Degree bin to sample from.
        bin_data (pd.DataFrame): DataFrame with gene degrees and bin assignments.
    Returns:
        int: Randomly selected gene ID from the specified bin.
    """
    possible_genes = bin_data[bin_data['bin'] == gene_bin].index
    return rn.choice(possible_genes)
        
def get_degree_bins(degrees):
    """Create degree bins with approximately 100 genes each for degree matching.
    Args:
        degrees (pd.DataFrame): DataFrame with node degrees (column 1 contains degree values).
    Returns:
        pd.DataFrame: Input DataFrame with added 'bin' column assigning each gene to a degree bin.
    """
    degree_counts = defaultdict(int)
    for x in degrees[1].values:
        degree_counts[x] += 1
    bins = [10000]
    bin_assignments = {}
    bin_totals = {}
    current_bin=10000
    current_total=0
    for degree_count in pd.DataFrame({'count':degree_counts}).index[::-1]:
        if current_total < 100:
            bin_assignments[degree_count] = current_bin
            current_total += degree_counts[degree_count]
        else:
            bin_totals[current_bin] = current_total
            current_bin = degree_count
            bins.append(current_bin)
            bin_assignments[degree_count] = current_bin
            current_total = degree_counts[degree_count]
    bin_totals[current_bin] = current_total
    degrees['bin'] = degrees[1].apply(lambda x: bin_assignments[x])
    return degrees

    
def write_simulated_geneset(geneset, outdir, setname, set_number, overlap, relevance, total_genes, repeat, background):
    """Write a simulated gene set to file with parameter-encoded filename.
    Args:
        geneset (list): List of gene IDs to write.
        outdir (str): Output directory for the file.
        setname (str): Base name for the gene set.
        set_number (str): Set identifier.
        overlap (int): Number of overlapping genes (for filename encoding).
        relevance (float): Relevance parameter (for filename encoding).
        total_genes (int): Total number of genes (for filename encoding).
        repeat (int): Repeat/iteration number (for filename encoding).
        background (str): Background matching method (for filename encoding).
    Returns:
        None. Writes gene list to file with newline-separated gene IDs.
    """
    outfile = os.path.join(outdir, f'{setname}_overlap{overlap}_relevance{relevance}_totalgenes{total_genes}_repeat{repeat}_background{background}_{set_number}.txt')
    with open(outfile, 'w') as out:
        out.write('\n'.join([str(x) for x in geneset]))
        
def write_simulated_geneset_quant(geneset, outdir, setname, set_number, overlap, relevance, total_genes, repeat, background):
    """Write a simulated gene set with scores to file with parameter-encoded filename.
    Args:
        geneset (dict): Dictionary mapping gene IDs to scores.
        outdir (str): Output directory for the file.
        setname (str): Base name for the gene set.
        set_number (str): Set identifier.
        overlap (int): Number of overlapping genes (for filename encoding).
        relevance (float): Relevance parameter (for filename encoding).
        total_genes (int): Total number of genes (for filename encoding).
        repeat (int): Repeat/iteration number (for filename encoding).
        background (str): Background matching method (for filename encoding).
    Returns:
        None. Writes gene profile to TSV file with 'Entrez' and 'P-value' columns.
    """
    outfile = os.path.join(outdir, f'{setname}_overlap{overlap}_relevance{relevance}_totalgenes{total_genes}_repeat{repeat}_background{background}_{set_number}.txt')
    out_df = pd.DataFrame({'P-value': geneset}).reset_index(names='Entrez')
    out_df.to_csv(outfile, sep='\t', index=False)
        

def check_genesets_against_network(genesets, network_node_file):
    """Filter gene sets to include only genes present in the network.
    Args:
        genesets (dict): Dictionary of gene sets (set name -> set of gene IDs).
        network_node_file (str): Path to file containing network node IDs.
    Returns:
        dict: Filtered gene sets containing only genes present in the network.
    """
    network_nodes = set(pd.read_csv(network_node_file, header=None, sep='\t')[0])
    for geneset in genesets:
        genesets[geneset] = genesets[geneset].intersection(network_nodes)
    return genesets



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--setfile', type=str, help='File containing all gene sets, or list of file names for gene set profiles')
    parser.add_argument('--outdir', type=str, help='')
    parser.add_argument('--netnodefile', type=str, help='Path to the file containing network nodes')
    parser.add_argument('--overlap', type=int, help='Number of overlapping genes')
    parser.add_argument('--relevance', type=float, help='Relevance of the gene set to control dilution')
    parser.add_argument('--totalgenes', type=int, help='total_genes')
    parser.add_argument('--nrepeats', type=int, help='How many interations to do for each parameter set')
    parser.add_argument('--background', type=str, help='')
    parser.add_argument('--quant', action='store_true', help='Should quantitative scores also be included? Note: only relevance implemented for quant.')
    args = parser.parse_args()

    if args.quant:
        with open(args.setfile, 'r') as f:
            file_list = [x.strip() for x in f.readlines()]
        for geneset_file in tqdm(file_list):
            gene_profile1 = load_full_gene_profile(os.path.join(os.path.dirname(args.setfile), geneset_file+'_CV.txt'), p_th=1e-5, gene_col='Entrez', score_col='P-value', return_dict=True)
            gene_profile2 = load_full_gene_profile(os.path.join(os.path.dirname(args.setfile), geneset_file+'_RV.txt'), p_th=1e-5, gene_col='Entrez', score_col='P-value', return_dict=True)
            for i in range(args.nrepeats):
                noised_set1 = add_noise_to_gene_set_quant(gene_profile1, args.relevance, args.background, list(gene_profile1.keys()), args.netnodefile)
                noised_set2 = add_noise_to_gene_set_quant(gene_profile2, args.relevance, args.background, list(gene_profile2.keys()), args.netnodefile)
                write_simulated_geneset_quant(noised_set1, args.outdir, geneset_file, 'CV', args.overlap, args.relevance, args.totalgenes, i+1, args.background)
                write_simulated_geneset_quant(noised_set2, args.outdir, geneset_file, 'RV', args.overlap, args.relevance, args.totalgenes, i+1, args.background)
    else:    
        all_genesets = load_node_sets(args.setfile)
        filtered_genesets = check_genesets_against_network(all_genesets, args.netnodefile)
        filtered_genesets = {k:v for k,v in filtered_genesets.items() if len(v) >= args.totalgenes}
        bakground_data = [] # should I make this have some degree of degree matching? Or GO matching? I.e. similar number of annotations. 
    # Note: genesets should already be prepocessed to those genes that are in the network
        for geneset in tqdm(filtered_genesets):
            all_set_genes = list(filtered_genesets[geneset])
            for i in range(args.nrepeats):
                # create the overlap partitions
                overlap_set, set1_unique, set2_unique = partition_gene_set(all_set_genes, args.totalgenes, args.overlap)
                #print(f'Overlap: {len(overlap_set)}, Set1: {len(set1_unique)}, Set2: {len(set2_unique)}')
                all_subset_genes = set1_unique + set2_unique + overlap_set
                if args.relevance < 1:
                # TODO add noise to the unique sets gene_set, relevance, background, all_subset_genes, netnodefile=None
                    noised_set1_unique, noised_set2_unique = add_noise_to_gene_set(set1_unique, args.relevance, args.background, all_subset_genes, args.netnodefile), add_noise_to_gene_set(set2_unique, args.relevance, args.background, all_subset_genes, args.netnodefile)
                    final_set1 = overlap_set + noised_set1_unique
                    final_set2 = overlap_set + noised_set2_unique
                else:
                    final_set1 = overlap_set + set1_unique
                    final_set2 = overlap_set + set2_unique
                write_simulated_geneset(final_set1, args.outdir, geneset, 'CV', args.overlap, args.relevance, args.totalgenes, i+1, args.background)
                write_simulated_geneset(final_set2, args.outdir, geneset, 'RV', args.overlap, args.relevance, args.totalgenes, i+1, args.background)

    #python create_sim_genesets.py --setfile /cellar/users/snwright/Data/NetColocTest/Reference/go.test --outdir /cellar/users/snwright/Data/NetColocTest/inputs/GO/ --netnodefile /cellar/users/snwright/Data/NetColocTest/Reference/pcnet2_0_nodelist.txt --overlap 0 --relevance 0.5 --totalgenes 100 --nrepeats 2 --background degree