# File containing general functions for working with genesets


import pandas as pd


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


def load_seed_genes(trait, rare_or_common, inputdir):
    '''From do_netcoloc'''
    assert rare_or_common in ['rare', 'common']
    rvc = {'rare':'RV', 'common':'CV'}[rare_or_common]
    try:
        seeds = pd.read_csv(os.path.join(inputdir, f"{trait}_{rvc}.txt"), sep='\t')
        if 'P-value' in seeds.columns:
            seeds = seeds.sort_values(by='P-value')['Entrez'].unique().tolist()
        else:
            seeds = seeds['Entrez'].tolist()
    except KeyError:
        seeds = pd.read_csv(os.path.join(inputdir, f"{trait}_{rvc}.txt"), sep='\t', header=None)[0].tolist()
    return seeds


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
