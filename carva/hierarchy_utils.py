import pandas as pd
import numpy as np
import networkx as nx
import os
import ndex2 as ndex
import sys
cwd = os.getcwd()
if cwd not in sys.path:
    sys.path.append(cwd)
from network_utils import *
from geneset_utils import *
import cdapsutil
from gprofiler import GProfiler
gp = GProfiler("MyToolName/0.1")


## Load subnetworks
def load_subnetwork(subnet_dir, trait_pair, net_suff, return_cx=False, name=''):
    """ Load a subnetwork from the specified directory and trait pair.
    Args:
        subnet_dir (str): Directory containing the subnetwork files.
        trait_pair (str): Trait pair identifier
        net_suff (str): Suffix for the network file (e.g., 'all', 'input').
        return_cx (bool, optional): If True, returns the network in CX format, otherwise returns the network as a networkx graph. Defaults to False.
        name (str, optional): Name to assign to the network. Defaults to ''.
    Returns:
        G (networkx.Graph or ndex2.NiceCXNetwork): The loaded subnetwork
    """
    net_file = os.path.join(subnet_dir, f'{trait_pair}_subnetwork_{net_suff}.tsv')
    net_df = pd.read_csv(net_file, sep='\t') # edgelist
    node_df = pd.read_csv(os.path.join(subnet_dir, f'{trait_pair}_subnetwork_{net_suff}_node_attributes.tsv'), sep='\t') # node attributes
    node_df = clean_nodes(node_df)

    # Construct a graph object from the edgelist
    edge_cols = [c for c in net_df.columns if c not in ('source','target')]
    if len(edge_cols) == 0:
        G =nx.Graph()
    else:
        G = nx.from_pandas_edgelist(
        net_df,
        source='source',
        target='target',
        edge_attr=edge_cols,
        create_using=nx.Graph()
        )
    # Ensure that all nodes are present (even if the have degree of zero)
    G.add_nodes_from(node_df['Entrez'])

    # Add node attributes
    node_attr_map = node_df.set_index('Entrez').to_dict('index')
    nx.set_node_attributes(G, node_attr_map)
    if return_cx:
        G_cx = ndex.create_nice_cx_from_networkx(G)
        G_cx.add_network_attribute('name', f'{name} {net_suff}')
        return G_cx
    else:
        return G


def clean_nodes(node_df):
    """ Ensure node attributes are in the correct format and fill missing values.
    Args:
        node_df (pandas.DataFrame): DataFrame containing node attributes.
    Returns:
        pandas.DataFrame: Cleaned DataFrame with appropriate types and filled missing values.
    """
    node_df['name'] = node_df.Entrez.astype(str)
    node_df = node_df.fillna({'symbol_R': '', 'symbol_C' :'', 'pval_R': 1, 'pval_C': 1})
    node_df['symbol_R'] = node_df.symbol_R.astype(str)
    node_df['symbol_C'] = node_df.symbol_C.astype(str)
    node_df['symbol'] = node_df.apply(lambda x: x.symbol_R if x.rare else x.symbol_C if x.common else '', axis=1)
    node_df['logp'] = node_df.apply(lambda x: -1* np.log10(min(x.pval_R, x.pval_C)+ 1e-250), axis=1)
    return node_df


def load_subnetwork_edges(trait_pair, subnet_dir, suffix='all'):
    """ Just load the edge list
    Args:
        trait_pair (str): Trait pair identifier
        subnet_dir (str): Directory containing the subnetwork files.
        suffix (str, optional): Suffix for the network file (e.g., 'all', 'input'). Defaults to 'all'.
    Returns:
        pandas.DataFrame: DataFrame containing the edge list.
    """
    all_df = pd.read_csv(os.path.join(subnet_dir, f'{trait_pair}_subnetwork_{suffix}.tsv'), sep='\t')
    return all_df


def load_subnetwork_node_info(trait_pair, subnet_dir, suffix='all'):
    """ Just load the node attributes
    Args:
        trait_pair (str): Trait pair identifier
        subnet_dir (str): Directory containing the subnetwork files.
        suffix (str, optional): Suffix for the network file (e.g., 'all', 'input'). Defaults to 'all'.
    Returns:
        pandas.DataFrame: DataFrame containing the node attributes.
    """
    all_nodes = pd.read_csv(os.path.join(subnet_dir, f'{trait_pair}_subnetwork_{suffix}_node_attributes.tsv'), sep='\t')
    all_nodes = clean_nodes(all_nodes)
    return all_nodes


## Create hierarchy ------------------------------------------

def create_hierarchy(G_cx, verbose=True, filter_nodes=None):
    """ Create a hierarchy from the given network in CX format.
    Args:

        G_cx (ndex2.NiceCXNetwork): The input network in CX format.
        verbose (bool, optional): If True, prints hierarchy details. Defaults to True.
        filter_nodes (list, optional): List of nodes to filter by coloc_gene attribute. If
            provided, only nodes with coloc_gene=1 will be included in the hierarchy. Defaults to None. If none, all input nodes and coloc genes will be included.
    Returns:
        G_hier (networkx.Graph): The hierarchy graph.
    """
    if filter_nodes is not None:
        name = G_cx.get_name()
        G = G_cx.to_networkx()
        coloc_nodes = [node for node in G.nodes(data=True) if G.nodes[node]['coloc_gene']==1]
        Gsub = G.subgraph(coloc_nodes)
        G_cx = ndex.create_nice_cx_from_networkx(G)
        G_cx.add_network_attribute('name', name)
    cd = cdapsutil.CommunityDetection()
    # Run the community detection algorithm hiDef
    CX_hier = cd.run_community_detection(G_cx, algorithm='hidefv1.1beta',arguments={'--maxres':'10'})
    if verbose:
        describe_cx_hierarchy(CX_hier)
    G_hier = CX_hier.to_networkx()
    return G_hier
    
def describe_cx_hierarchy(G_hier):
    print('Hierarchy name: ' + str(G_hier.get_name()))
    print('# nodes: ' + str(len(G_hier.get_nodes())))
    print('# edges: ' + str(len(G_hier.get_edges())))


## Annotate hierarchy --------------------------------------
def create_hier_df(G_hier):
    """ Create a DataFrame from the hierarchy graph.
    Args:
        G_hier (networkx.Graph): The hierarchy graph.
    Returns:
        pandas.DataFrame: DataFrame containing hierarchy node attributes.
    """
    hier_df = pd.DataFrame.from_dict(dict(G_hier.nodes(data=True)), orient="index")
    hier_nodes = {}
    for node in G_hier.nodes(data=True):
        gene_list = [int(x) for x in node[1]['CD_MemberList'].split(' ')]
        symbol_list = map_genes_using_network_map(gene_list)
        hier_nodes[node[0]] = ' '.join(list(symbol_list.values()))
    hier_df['SymbolList'] = hier_df.index.map(hier_nodes)
    for col in ['CD_MemberList_Size', 'CD_MemberList_LogSize', 'HiDeF_persistence']:
        hier_df[col] = hier_df[col].astype(float)
    return hier_df.loc[:, ('CD_MemberList', 'SymbolList', 'CD_MemberList_Size', 'CD_MemberList_LogSize', 'HiDeF_persistence')]


def name_hierarchy_systems(hier_df, outdir=None, gene_col='CD_MemberList', write=False, hier_name=None, sources=['GO:BP'],
                           min_term_size=50, max_term_size=1000, p_value_threshold=1e-3, min_intersection=3):
    """ Annotate the hierarchy with GO terms for each community.
    Args:
        hier_df (pandas.DataFrame): DataFrame containing hierarchy node attributes.
        outdir (str, optional): Directory to write the results to. Defaults to None.
        gene_col (str, optional): Column name containing the gene lists for each community. Defaults to 'CD_MemberList'.
        write (bool, optional): If True, writes the results to a file. Defaults to False.
        hier_name (str, optional): Name for the hierarchy, used in output file names. Defaults to None.
        sources (list, optional): List of GO sources to use for annotation. Defaults to ['GO:BP'].
        min_term_size (int, optional): Minimum size of a GO term to keep. Defaults to 50.
        max_term_size (int, optional): Maximum size of a GO term to keep. Defaults to 1000.
        p_value_threshold (float, optional): Significance threshold for GO terms. Defaults to 1e-3.
        min_intersection (int, optional): Minimum number of community terms annotated to the GO term. Defaults to 3.
    Returns:
        pandas.DataFrame: Annotated hierarchy DataFrame with GO names.
    """
    system_names = {}
    all_results = []
    for comm in hier_df.index.tolist():
        # get the genes in the community
        if isinstance(hier_df[gene_col].loc[comm], str):
            focal_genes = hier_df[gene_col].loc[comm].split(' ')
        else:
            focal_genes = hier_df[gene_col].loc[comm]
        system_names[comm] = comm # default is to keep identifier as name, will be replaced if confident annotation available
        if len(focal_genes) > 2:
            # get all available GO annotations
            gp_temp = pd.DataFrame(gp.profile(focal_genes,significance_threshold_method='fdr',
                                                sources=sources, no_evidences=False))
            if len(gp_temp) > 0: # make sure data is not empty
                # filter and append the results
                gp_temp["community"] = comm
                all_results.append(gp_temp)
                gp_temp = filter_go_annotations(gp_temp, min_term_size, max_term_size, p_value_threshold, min_intersection)
                if len(gp_temp) > 1:
                    system_names[comm] = gp_temp.head(1)['name'].tolist()[0]
    go_results = pd.concat(all_results)
    if write:
        assert outdir is not None, "Please provide an output directory to write the results to."
        assert hier_name is not None, "Please provide a name for the hierarchy."
        go_results.to_csv(os.path.join(outdir, hier_name+"_go_results.tsv"), sep="\t")
    hier_df['GO_Name'] = hier_df.index.map(system_names)
    return hier_df

def filter_go_annotations(go_df, term_min=50, term_max=1000, p_th=1e-4, min_intersection=3):
    """Filters available annotations for a community based on specificity and significance.
    Args:
        go_df (pandas.DataFrame): All available significant GO terms for each community
        term_min (int, optional):   The minimum size of a term to keep. Defaults to 50.
        term_max (int, optional): The maximum size of a term to keep. Defaults to 1000.
        p_th (float, optional): The significance threshold. Defaults to 1e-4.
        min_intersection (int, optional): Minimum number of community terms annotated to the GO term. Defaults to 3.

    Returns:
        pandas.DataFrame: A filter dataframe of GO annotations per community, sorted by sum of precision and recall.
    """
    go_df = go_df[(go_df['term_size'] <= term_max) & (go_df['term_size'] >= term_min)]
    go_df = go_df[go_df['intersection_size'] >= min_intersection]
    go_df = go_df[go_df['p_value'] < p_th] # set a stringent pvalue threshold
    go_df['sum_PR'] = go_df['recall'] + go_df['precision']
    go_df = go_df.sort_values('sum_PR',ascending=False)
    return go_df


def add_seed_gene_fractions(hier_df, G_cx):
    """ Calculate the fraction of seed genes in each community and their average z-scores.

    Args:
        hier_df (pandas.DataFrame): DataFrame containing hierarchy node attributes.
        G_cx (ndex2.NiceCXNetwork): The input network in CX format
    Returns:
        pandas.DataFrame: Updated hierarchy DataFrame with seed gene fractions and z-scores.
    """
    seed_fracs = {}
    for comm in hier_df.index.tolist():
        gene_list = [int(x) for x in hier_df['CD_MemberList'].loc[comm].split(' ')]
        rare_count = 0
        common_count = 0
        shared_count = 0
        rare_z = 0
        common_z = 0
        shared_z = 0
        n = len(gene_list)
        for gene in gene_list:
            if G_cx.get_node_attribute(gene,'rare')['v'] == 'true':
                rare_count += 1
            if G_cx.get_node_attribute(gene,'common')['v'] == 'true':
                common_count += 1
            if G_cx.get_node_attribute(gene,'shared')['v'] == 'true':
                shared_count += 1
            rare_z += float(G_cx.get_node_attribute(gene,'z_R')['v'])
            common_z += float(G_cx.get_node_attribute(gene,'z_C')['v'])
            shared_z += float(G_cx.get_node_attribute(gene,'Z_coloc')['v'])
        seed_fracs[comm] = {'rare': rare_count/n, 'common': common_count/n, 'shared': shared_count/n,
                            'rare_z': rare_z/n, 'common_z': common_z/n, 'shared_z': shared_z/n, 'rc_ratio': -1 * np.log2(rare_z/common_z)}
    seed_df = pd.DataFrame.from_dict(seed_fracs, orient='index')
    hier_df = hier_df.join(seed_df, how='inner')
    return hier_df
    
### I/O functions ----------------------------------------

def write_nx_hierarchy(G_hier, hier_df, outdir,name):
    """ Write the hierarchy graph and its attributes to files.
    Args:
        G_hier (networkx.Graph): The hierarchy graph.
        hier_df (pandas.DataFrame): DataFrame containing hierarchy node attributes.
        outdir (str): Directory to write the results to.
        name (str): Name for the hierarchy, used in output file names.
    """
    hier_df.to_csv(os.path.join(outdir, name+"_hierarchy_nodes.tsv"), sep="\t")
    nx.to_pandas_edgelist(G_hier).to_csv(os.path.join(outdir, name+"_hierarchy_edges.tsv"), sep="\t", index=False)

def add_annotations_to_hierarchy(G_hier, hier_df, annot_cols=['GO_Name', 'rare', 'common', 'shared', 'rare_z', 'common_z', 'shared_z', 'rc_ratio','CD_MemberList',
                                                              'SymbolList', 'CD_MemberList_Size', 'CD_MemberList_LogSize', 'HiDeF_persistence']):
    """ Add annotations to the hierarchy graph from the hierarchy DataFrame.
    Args:
        G_hier (networkx.Graph): The hierarchy graph.
        hier_df (pandas.DataFrame): DataFrame containing hierarchy node attributes.
        annot_cols (list): List of columns to annotate in the graph.
    """
    G_hier_out = nx.Graph()
    G_hier_out.add_nodes_from(G_hier.nodes())
    G_hier_out.add_edges_from(G_hier.edges())
    for col in annot_cols:
        nx.set_node_attributes(G_hier_out, hier_df[col].to_dict(), name=col)
    return G_hier_out

def upload_cx_hierarchy(G_hier, hier_df, name, username, password, template_uuid=None,
                        annot_cols=['GO_Name', 'rare', 'common', 'shared', 'rare_z', 'common_z', 'shared_z', 'rc_ratio','CD_MemberList', 
                                    'SymbolList', 'CD_MemberList_Size','CD_MemberList_LogSize', 'HiDeF_persistence']):
    """ Add annotations and upload the hierarchy graph to NDEx.
    Args:
        G_hier (networkx.Graph): The hierarchy graph.
        hier_df (pandas.DataFrame): DataFrame containing hierarchy node attributes.
        name (str): Name for the hierarchy, used in NDEx network name.
        username (str): NDEx username.
        password (str): NDEx password.
        template_uuid (str, optional): UUID of the template to use for the network. Defaults to None.
        annot_cols (list, optional): List of columns to annotate in the graph. Defaults to ['GO_Name', 'rare', 'common', 'shared', 'rare_z', 'common_z', 'shared_z', 'rc_ratio','CD_MemberList',
                                    'SymbolList', 'CD_MemberList_Size','CD_MemberList_LogSize', 'HiDeF_persistence'].
    Returns:
        str: NDEx network UUID of the uploaded hierarchy.
    """
    G_out = add_annotations_to_hierarchy(G_hier, hier_df, annot_cols=annot_cols)
    # apply style
    id_str = upload_network(G_out, name+'_Hierarchy', username=username, password=password, template=template_uuid)
    return id_str
