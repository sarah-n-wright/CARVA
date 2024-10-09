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

def partition_gene_set(all_genes, total_genes, overlap):
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

        
def get_matched_gene(gene_bin, bin_data):
    possible_genes = bin_data[bin_data['bin'] == gene_bin].index
    return rn.choice(possible_genes)
        
def get_degree_bins(degrees):
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
    outfile = os.path.join(outdir, f'{setname}_overlap{overlap}_relevance{relevance}_totalgenes{total_genes}_repeat{repeat}_background{background}_{set_number}.txt')
    with open(outfile, 'w') as out:
        out.write('\n'.join([str(x) for x in geneset]))
        

def load_node_sets(node_set_file, delimiter='\t', verbose=False, id_type="Entrez"):
    """ Load node sets from a text file into a dictionary
    
    Args:
        node_set_file (str): path to node set file
        delimiter (str): delimiter for node set file
        verbose (bool): print out number of node sets loaded
        id_type (str): type of node ID to use for graph
    
    Returns:
        dict: dictionary of node sets
    """
    f = open(node_set_file)
    node_set_lines = f.read().splitlines()
    node_set_lines_split = [line.split(delimiter) for line in node_set_lines]
    f.close()
    node_sets = {node_set[0]:set(node_set[1:]) for node_set in node_set_lines_split}
    if id_type == "Entrez":
        for set_id in node_sets:
            node_sets[set_id] = {int(node) for node in list(node_sets[set_id]) if node.isnumeric()}
    if verbose:
        print('Node cohorts loaded:', node_set_file)
    return node_sets

def check_genesets_against_network(genesets, network_node_file):
    network_nodes = set(pd.read_csv(network_node_file, header=None, sep='\t')[0])
    for geneset in genesets:
        genesets[geneset] = genesets[geneset].intersection(network_nodes)
    return genesets



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--setfile', type=str, help='File containing all gene sets')
    parser.add_argument('--outdir', type=str, help='')
    parser.add_argument('--netnodefile', type=str, help='Path with the reference gene sets')
    parser.add_argument('--overlap', type=int, help='Number of overlapping genes')
    parser.add_argument('--relevance', type=float, help='Relevance of the gene set to control dilution')
    parser.add_argument('--totalgenes', type=int, help='total_genes')
    parser.add_argument('--nrepeats', type=int, help='How many interations to do for each parameter set')
    parser.add_argument('--background', type=str, help='')
    args = parser.parse_args()
    
    
    
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