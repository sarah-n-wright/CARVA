# File containing general functions for working with genesets
import pandas as pd
import os
import numpy as np
CWD = os.path.dirname(os.path.abspath(__file__))

# TODO reduce redundancy of gene set loading functions. 

def map_genes_using_network_map(genelist, network_map_file=os.path.join(CWD, '../Reference_Data/pcnet2_0_node_map.txt'),
                               from_col='Entrez', to_col='Symbol'):
    """ Map genes from one ID type to another using a network map file.
    Args:
        genelist (list): List of gene IDs to map.
        network_map_file (str): Path to the network map file.
        from_col (str): Column name in the map file for the original gene IDs.
        to_col (str): Column name in the map file for the target gene IDs.
    Returns:
        dict: Dictionary mapping original gene IDs to target gene IDs.
    """
    map_df = pd.read_csv(network_map_file, sep='\t')
    genemap = map_df.set_index(from_col)[to_col].to_dict()
    
    return {g: genemap[g] for g in genelist if g in genemap}


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


def load_seed_genes(trait, rare_or_common, inputdir, usecol='Entrez'):
    """ Load specified gene set
    Args:
        trait (str): Trait name to load seeds for.
        rare_or_common (str): 'rare' or 'common' to specify the type of seed genes.
        inputdir (str): Directory containing the seed files.
        usecol (str): Column to use for loading seeds, defaults to 'Entrez'.
    Returns:
        list: List of unique seed genes.
    """
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
    """ Load a gene profile from a file, filtering by p-value threshold and minimum number of genes.
    Args:
        file (str): Path to the gene profile file.
        p_th (float): P-value threshold for filtering genes.
        min_genes (int): Minimum number of genes required to return a list.
    Returns:
        list: List of unique gene IDs that meet the criteria, or None if not enough genes
        """
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
        id_type (str): Type of ID to use for output, e.g., 'Entrez' or 'Symbol'

    Usage:
        my_dict = {'key1': [1, 2, 3], 'key2': ['a', 'b', 'c']}
        write_dict_to_file(my_dict, 'output.tsv')
    """
    with open(output_file, 'w') as f:
        for key, values in gene_dict.items():
            # Convert all values to strings and join with tabs
            if id_type != 'Entrez':
                line = f"{key}\t" + "\t".join(str(v) for v in values)
            else:
                line = f"{key}\t" + "\t".join(str(int(v)) for v in values)
            f.write(line + "\n")
    
    
def load_full_gene_profile(file, p_th, gene_col='Entrez', score_col='P-value', return_dict=True):
    """ Load a full gene profile from a file, filtering by p-value threshold.
    Args:
        file (str): Path to the gene profile file.
        p_th (float): P-value threshold for filtering genes.
        gene_col (str): Column name for gene IDs.
        score_col (str): Column name for gene scores.
        return_dict (bool): Whether to return a dictionary or a DataFrame.
    Returns:
        dict or pandas.DataFrame: If return_dict is True, returns a dictionary mapping gene IDs to scores.
                                  If False, returns a DataFrame with gene IDs and scores.   
    """
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
    """ Partition a gene set into two equal-size subsets.
    Args:
        filepath (str): Path to the gene profile file.
        outpref (str): Output prefix for the split files.
        p_th (float, optional): P-value threshold for filtering genes. 
        min_genes (int): Minimum number of genes required in each subset.
        gene_col (str): Column name for gene IDs in the input file.
        score_col (str): Column name for gene scores in the input file.
        write_results (bool): Whether to write the results to files or return them as DataFrames.
    Returns:
        None or tuple: If write_results is True, writes the results to files and returns None.
                       If False, returns a tuple of two DataFrames containing the split gene sets.
    """
    if (p_th is not None) and (p_th < 0.5):
        prof = load_full_gene_profile(filepath, p_th, gene_col=gene_col, score_col=score_col, return_dict=False)
        n_genes = len(prof)
        if n_genes < min_genes:
            print(f'Not enough genes at p_th {p_th} for profile {filepath}.')
            return None
        else:
        # split the sig genes in half randomly
            genes1 = prof.sample(frac=0.5)
            genes2 = prof.drop(genes1.index)
            if write_results:
                genes1.to_csv(outpref + os.path.basename(filepath).split('.txt')[0] + '_RV.txt', sep='\t', index=False)
                genes2.to_csv(outpref + os.path.basename(filepath).split('.txt')[0] + '_CV.txt', sep='\t', index=False)
            else:
                return genes1, genes2
    elif (p_th is None) or (p_th > 0.5):
        # randomize half the genes using uniform distribution of p-values which is the null hypothesis. 
        prof = load_full_gene_profile(filepath, p_th = 1, gene_col=gene_col, score_col=score_col, return_dict=False)
        # create rare_gene profile
        genes1 = prof.sample(frac=0.5)
        genes1_random = prof.drop(genes1.index)
        genes1_random[score_col] = np.random.uniform(0, 1, size=len(genes1_random))
        genes1_out = pd.concat([genes1, genes1_random])

        # create common_gene profile
        genes2 = prof.drop(genes1.index)
        genes2_random = genes1.copy()
        genes2_random[score_col] = np.random.uniform(0, 1, size=len(genes2_random))
        genes2_out = pd.concat([genes2, genes2_random])

        if write_results:
            genes1_out.to_csv(outpref+ os.path.basename(filepath).split('.txt')[0] + '_RV.txt', sep='\t', index=False)
            genes2_out.to_csv(outpref + os.path.basename(filepath).split('.txt')[0] + '_CV.txt', sep='\t', index=False)
        else:
            return genes1_out, genes2_out
