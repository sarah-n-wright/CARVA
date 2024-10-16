# File containing general functions for workign with networks

import sys
import pandas as pd
import os
import numpy as np
import networkx as nx
import ndex2
from getpass import getpass

def load_saved_network_nodes(indir, net_name):
    '''From do_netcoloc'''
    if exists(os.path.join(indir, net_name+ "_nodes.txt")):
        pc_nodes = pd.read_csv(os.path.join(indir, net_name+ "_nodes.txt"), sep='\t', header=None, index_col=0).index.tolist()
        pc_nodes = [int(x) for x in pc_nodes]
        return pc_nodes
    else:
        return
        
def load_saved_network_degrees(indir, net_name):
    '''From do_netcoloc'''
    if exists(os.path.join(indir, net_name+ "_degrees.txt")):
        degree_map = pd.read_csv(os.path.join(indir, net_name+ "_degrees.txt"), sep='\t', header=None, index_col=0).to_dict()[1]
        return degree_map
    else:
        return
        
def create_saved_nodes_and_degrees(uuid, outdir, net_name, nodes=True, degrees=True):
    '''From do_netcoloc'''
    G_PC = load_network(uuid)
    if nodes:
        pc_nodes_df = pd.DataFrame({'node': G_PC.nodes})
        pc_nodes_df['Entrez'] = pc_nodes_df.node.apply(lambda x: G_PC.nodes[x]['GeneID'])
        pd.DataFrame(pc_nodes_df.Entrez).to_csv(os.path.join(outdir, net_name+ "_nodes.txt"), sep='\t', header=False, index=False)
    if degrees:
        degree_map = pd.DataFrame(G_PC.degree())
        degree_map['Entrez'] = degree_map[0].apply(lambda x: G_PC.nodes[x]['GeneID'])
        degree_map.loc[:, ('Entrez', 1)].to_csv(os.path.join(outdir, net_name+ "_degrees.txt"), sep='\t', header=False, index=False)


def load_network(uuid='e8cc9239-d91a-11eb-b666-0ac135e8bacf', use_password=False, verbose=True):
    """Wrapper function for loading a network from NDEx

    Args:
        uuid (str, optional): NDEx identifier. Defaults to 'e8cc9239-d91a-11eb-b666-0ac135e8bacf'.
        use_password (bool, optional): set True if the network is not public. Defaults to False.

    Returns:
        :py:class:`networkx.Graph`: A networkx object of the desired network
    """
    ndex_server='public.ndexbio.org'
    if use_password:
        ndex_user=getpass.getpass("NDEX Username:")
        ndex_password=getpass.getpass("NDEX Password:")
    else:
        ndex_user=None
        ndex_password=None
    G_cx = ndex2.create_nice_cx_from_server(
            ndex_server, 
            username=ndex_user, 
            password=ndex_password, 
            uuid=uuid)
    G = G_cx.to_networkx()
    if verbose:
        print(f'Network Name:{G_cx.get_name()}')
        print(f'Number of nodes: {len(G_overlap.nodes)}')
        print(f'Number of edges: {len(G_overlap.edges)}')
    return G