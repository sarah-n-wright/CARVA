# File containing general functions for working with genesets


import pandas as pd
import os
import numpy as np


def map_genes_using_network_map(genelist, network_map_file='/cellar/users/snwright/Data/RareCommon/inputs/pcnet2_0_node_map.txt',
                               from_col='Entrez', to_col='Symbol'):
    map_df = pd.read_csv(network_map_file, sep='\t')
    genemap = map_df.set_index(from_col)[to_col].to_dict()
    
    return {g: genemap[g] for g in genelist if g in genemap}


def load_node_sets(node_set_file, delimiter='\t', verbose=False, id_type="Entrez"):
    """ Load node sets from a text file into a dictionary
    [From create_sim_genesets] 
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


def load_seed_genes(trait, rare_or_common, inputdir, usecol='Entrez'):
    '''From do_netcoloc'''
    assert rare_or_common in ['rare', 'common']
    rvc = {'rare':'RV', 'common':'CV'}[rare_or_common]
    seeds = pd.read_csv(os.path.join(inputdir, f"{trait}_{rvc}.txt"), sep='\t').drop_duplicates()
    if usecol == 'Entrez':
        if 'P-value' in seeds.columns:
            seeds = seeds.sort_values(by='P-value')['Entrez'].unique().tolist()
        else:
            seeds = seeds['Entrez'].tolist()
    elif usecol in seeds.columns:
        seeds = seeds[usecol].tolist()
    else:
        seeds = pd.read_csv(os.path.join(inputdir, f"{trait}_{rvc}.txt"), sep='\t', header=None)[0].tolist()
    return list(set(seeds))


def load_gene_profile(file, p_th, min_genes):
    '''From gene_overlap'''
    gene_df = pd.read_csv(file, sep='\t')
    if 'P-value' in gene_df.columns:
        gene_df = gene_df[gene_df['P-value'] < p_th]
        gene_list = gene_df['Entrez'].unique().tolist()
    else:
        gene_df = pd.read_csv(file, sep='\t', header=None)
        gene_list = gene_df[0].tolist()
    if len(gene_list) < min_genes:
        return None
    else:
        return gene_list
    
    
def write_node_sets(gene_dict, output_file, id_type='Symbol'):
    """
    Write a dictionary to a file where each line represents one entry with the key first,
    followed by tab-separated list values.
    
    Args:
        dictionary (dict): Dictionary to write to file
        output_file (str): Path to output file
    """
    with open(output_file, 'w') as f:
        for key, values in gene_dict.items():
            # Convert all values to strings and join with tabs
            if id_type != 'Entrez':
                line = f"{key}\t" + "\t".join(str(v) for v in values)
            else:
                line = f"{key}\t" + "\t".join(str(int(v)) for v in values)
            f.write(line + "\n")

# Example usage:
# my_dict = {'key1': [1, 2, 3], 'key2': ['a', 'b', 'c']}
# write_dict_to_file(my_dict, 'output.tsv')
    
    
def load_full_gene_profile(file, p_th, gene_col='Entrez', score_col='P-value', return_dict=True):
    gene_df = pd.read_csv(file, sep='\t')
    assert gene_col in gene_df.columns, f"{gene_col} not in {file}"
    assert score_col in gene_df.columns, f"{score_col} not in {file}"
    gene_df = gene_df[gene_df[score_col] <= p_th].loc[:, [gene_col, score_col]].drop_duplicates()
    # create score dictionary
    if return_dict:
        gene_dict = {gene:score for gene, score in zip(gene_df[gene_col], gene_df[score_col])}
        return gene_dict
    else:
        return gene_df.loc[:, [gene_col, score_col]]

def split_gene_profile(filepath, outpref, p_th=None, min_genes=6, gene_col='Entrez', score_col='P-value', write_results=True):
    if (p_th is not None) and (p_th < 0.5):
        prof = load_full_gene_profile(filepath, p_th, gene_col=gene_col, score_col=score_col, return_dict=False)
        n_genes = len(prof)
        if n_genes < min_genes:
            print(f'Not enough genes at p_th {p_th} for profile {filepath}.')
            return None
        else:
        # split the sig genes in half randomly
            rare_genes = prof.sample(frac=0.5)
            common_genes = prof.drop(rare_genes.index)
            if write_results:
                rare_genes.to_csv(outpref + os.path.basename(filepath).split('.txt')[0] + '_RV.txt', sep='\t', index=False)
                common_genes.to_csv(outpref + os.path.basename(filepath).split('.txt')[0] + '_CV.txt', sep='\t', index=False)
            else:
                return rare_genes, common_genes
    elif (p_th is None) or (p_th > 0.5):
        # randomize half the genes using uniform distribution of p-values which is the null hypothesis. 
        prof = load_full_gene_profile(filepath, p_th = 1, gene_col=gene_col, score_col=score_col, return_dict=False)
        # create rare_gene profile
        rare_genes = prof.sample(frac=0.5)
        rare_random = prof.drop(rare_genes.index)
        rare_random[score_col] = np.random.uniform(0, 1, size=len(rare_random))
        rare_genes_out = pd.concat([rare_genes, rare_random])

        # create common_gene profile
        common_genes = prof.drop(rare_genes.index)
        common_random = rare_genes.copy()
        common_random[score_col] = np.random.uniform(0, 1, size=len(common_random))
        common_genes_out = pd.concat([common_genes, common_random])
        
        if write_results:
            rare_genes_out.to_csv(outpref+ os.path.basename(filepath).split('.txt')[0] + '_RV.txt', sep='\t', index=False)
            common_genes_out.to_csv(outpref + os.path.basename(filepath).split('.txt')[0] + '_CV.txt', sep='\t', index=False)
        else:
            return rare_genes_out, common_genes_out    
    