import pandas as pd
import numpy as np
import networkx as nx
import os
import sys
sys.path.append('/cellar/users/snwright/Git/rare_common/carva')
from network_utils import *
from geneset_utils import *
from getpass import getpass
import cdapsutil
from gprofiler import GProfiler
gp = GProfiler("MyToolName/0.1")

## Create hierarchy ------------------------------------------

def create_hierarchy(G_cx, verbose=True):
    cd = cdapsutil.CommunityDetection()
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


def name_hierarchy_systems(hier_df, outdir=None, gene_col='CD_MemberList', write=False, hier_name=None):
    system_names = {}
    all_results = []
    for comm in hier_df.index.tolist():
        # get the genes in the community
        if isinstance(hier_df[gene_col].loc[comm], str):
            focal_genes=hier_df[gene_col].loc[comm].split(' ')
        else:
            focal_genes=hier_df[gene_col].loc[comm]
        #print(comm,":", len(focal_genes))
        system_names[comm] = comm # default is to keep identifier as name, will be replaced if confident annotation available
        if len(focal_genes)>2:
            # get all available GO:BP annotations
            gp_temp = pd.DataFrame(gp.profile(focal_genes,significance_threshold_method='fdr',
                                                sources=['GO:BP'], no_evidences=False))
            if len(gp_temp)>0: # make sure data is not empty
                # filter and append the results
                gp_temp["community"] = comm
                all_results.append(gp_temp)
                gp_temp = filter_go_annotations(gp_temp, 50, 1000, 1e-3, 3)
                if len(gp_temp)>1:
                    system_names[comm] = gp_temp.head(1)['name'].tolist()[0]
                    #all_results.append(gp_temp)
    go_results = pd.concat(all_results)
    if write:
        assert outdir is not None, "Please provide an output directory to write the results to."
        assert hier_name is not None, "Please provide a name for the hierarchy."
        go_results.to_csv(os.path.join(outdir, hier_name+"_go_results.tsv"), sep="\t")
    #go_results.to_csv(outdir+"test.csv", index=True, sep="\t")
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

def filter_hierachy():
    pass

def add_seed_gene_fractions(hier_df, G_cx):
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
    hier_df.to_csv(os.path.join(outdir, name+"_hierarchy_nodes.tsv"), sep="\t")
    nx.to_pandas_edgelist(G_hier).to_csv(os.path.join(outdir, name+"_hierarchy_edges.tsv"), sep="\t", index=False)

def add_annotations_to_hierarchy(G_hier, hier_df, annot_cols=['GO_Name', 'rare', 'common', 'shared', 'rare_z', 'common_z', 'shared_z', 'rc_ratio','CD_MemberList',
                                                              'SymbolList', 'CD_MemberList_Size', 'CD_MemberList_LogSize', 'HiDeF_persistence']):
    G_hier_out = nx.Graph()
    G_hier_out.add_nodes_from(G_hier.nodes())
    G_hier_out.add_edges_from(G_hier.edges())
    for col in annot_cols:
        nx.set_node_attributes(G_hier_out, hier_df[col].to_dict(), name=col)
    return G_hier_out

def upload_cx_hierarchy(G_hier, hier_df, template_uuid, name, username, password, annot_cols=['GO_Name', 'rare', 'common', 'shared', 'rare_z', 'common_z', 'shared_z', 'rc_ratio','CD_MemberList', 'SymbolList', 'CD_MemberList_Size','CD_MemberList_LogSize', 'HiDeF_persistence']):
    G_out = add_annotations_to_hierarchy(G_hier, hier_df, annot_cols=annot_cols)

    # apply layout

    # apply style
    id_str = upload_network(G_out, name+'_Hierarchy', username=username, password=password, template=template_uuid)
    return id_str
