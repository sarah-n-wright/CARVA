# File containing general functions for workign with networks

import sys
import pandas as pd
import os
import numpy as np
import networkx as nx
import ndex2
from getpass import getpass

def upload_network(G, name, template=None, networkset=None, networkset_name=None, username=None, password=None):
    if (username is None) or (password is None):
        username=getpass('Username:')
        password = getpass('Password:')
    client = ndex2.client.Ndex2(username=username, password=password)

    G_cx = ndex2.create_nice_cx_from_networkx(G)
    G_cx.add_network_attribute('name', name)
    
    SERVER = 'http://ndexbio.org'
    if template is not None:
        G_cx.apply_template(SERVER, template, username=username, password=password)
    id_str = G_cx.upload_to(SERVER, username=username, password=password)
    id_str = id_str.split('/')[-1]
    if networkset is None:
        return id_str
    elif networkset =='create':
        set_id = client.create_networkset(networkset_name, description=networkset_name)
        set_id = set_id.split('/')[-1]
        client.add_networks_to_networkset(set_id, [id_str])
        return (set_id, id_str)
    else:
        client.add_networks_to_networkset(networkset, [id_str])
        return id_str

    
    
def load_saved_network_nodes(indir, net_name):
    '''From do_netcoloc'''
    if os.path.exists(os.path.join(indir, net_name+ "_nodes.txt")):
        pc_nodes = pd.read_csv(os.path.join(indir, net_name+ "_nodes.txt"), sep='\t', header=None, index_col=0).index.tolist()
        pc_nodes = [int(x) for x in pc_nodes]
        return pc_nodes
    else:
        return

def load_saved_network_degrees(indir, net_name):
    '''From do_netcoloc'''
    if os.path.exists(os.path.join(indir, net_name+ "_degrees.txt")):
        degree_map = pd.read_csv(os.path.join(indir, net_name+ "_degrees.txt"), sep='\t', header=None, index_col=0).to_dict()[1]
        return degree_map
    else:
        return

def create_saved_nodes_and_degrees(uuid, outdir, net_name, nodes=True, degrees=True, represents=False, alias='ncbigene:'):
    '''From do_netcoloc'''
    G_PC = load_network(uuid)
    if nodes:
        pc_nodes_df = pd.DataFrame({'node': G_PC.nodes})
        if not represents:
            pc_nodes_df['Entrez'] = pc_nodes_df.node.apply(lambda x: G_PC.nodes[x]['GeneID'])
        elif alias is not None:
            pc_nodes_df['Entrez'] = pc_nodes_df.node.apply(lambda x: [z for z in G_PC.nodes[x]['alias'] if alias in z][0].split(alias)[0])
        else:
            pc_nodes_df['Entrez'] = pc_nodes_df.node.apply(lambda x: G_PC.nodes[x]['represents'].split('ncbigene:')[-1])
        pd.DataFrame(pc_nodes_df.Entrez).to_csv(os.path.join(outdir, net_name+ "_nodes.txt"), sep='\t', header=False, index=False)
    if degrees:
        degree_map = pd.DataFrame(G_PC.degree())
        if not represents:
            degree_map['Entrez'] = degree_map[0].apply(lambda x: G_PC.nodes[x]['GeneID'])
        else:
            degree_map['Entrez'] = degree_map[0].apply(lambda x: G_PC.nodes[x]['represents'].split('ncbigene:')[-1])
        degree_map.loc[:, ('Entrez', 1)].to_csv(os.path.join(outdir, net_name+ "_degrees.txt"), sep='\t', header=False, index=False)


def create_ncbi_symbol_map_from_network(uuid, outdir):
    G = load_network(uuid)
    node_map = {x[1]['GeneID']: x[0] for x in G.nodes(data=True)}
    node_df = pd.DataFrame({'Symbol':node_map}).reset_index().rename(columns={'index':'Entrez'})
    node_df.to_csv(os.path.join(outdir, f'{net_name}_node_map.txt'), sep='\t', index=False)  
        
# +
import sys
import pandas as pd
import os
import numpy as np
import networkx as nx
import ndex2
from getpass import getpass

def load_network(uuid='e8cc9239-d91a-11eb-b666-0ac135e8bacf', use_password=False, verbose=True, 
                return_cx=False, ndex_password=None, ndex_user=None):
    """Wrapper function for loading a network from NDEx
s
    Args:
        uuid (str, optional): NDEx identifier. Defaults to 'e8cc9239-d91a-11eb-b666-0ac135e8bacf'.
        use_password (bool, optional): set True if the network is not public to be prompted to enter username and password.
                                        Default is False for public networks. 

    Returns:
        :py:class:`networkx.Graph`: A networkx object of the desired network
    """
    ndex_server='public.ndexbio.org'
    if use_password:
        if ndex_user is None:
            ndex_user=getpass("NDEX Username:")
        if ndex_password is None:
            ndex_password=getpass("NDEX Password:")
    else:
        ndex_user=None
        ndex_password=None
    G_cx = ndex2.create_nice_cx_from_server(
            ndex_server, 
            username=ndex_user, 
            password=ndex_password, 
            uuid=uuid)
    if return_cx:
        return G_cx
    G = G_cx.to_networkx()
    if verbose:
        print(f'Network Name:{G_cx.get_name()}')
        print(f'Number of nodes: {len(G.nodes)}')
        print(f'Number of edges: {len(G.edges)}')
    return G
# -