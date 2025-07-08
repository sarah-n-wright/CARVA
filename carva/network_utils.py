# File containing general functions for working with networks
import pandas as pd
import os
import ndex2
from getpass import getpass

def load_network(uuid='e8cc9239-d91a-11eb-b666-0ac135e8bacf', use_password=False, verbose=True, 
                return_cx=False, ndex_password=None, ndex_user=None):
    """Wrapper function for loading a network from NDEx
    Args:
        uuid (str, optional): NDEx identifier. Defaults to 'e8cc9239-d91a-11eb-b666-0ac135e8bacf' for PCNET
        use_password (bool, optional): set True if the network is not public to be prompted to enter username and password.
                                        Default is False for public networks. 
        verbose (bool, optional): If True, prints network name and number of nodes and edges. Defaults to True.
        return_cx (bool, optional): If True, returns the CX object instead of a networkx object. Defaults to False.
        ndex_password (str, optional): Password for NDEx. If not provided and use_password is True, will prompt for password.
        ndex_user (str, optional): Username for NDEx. If not provided and use_password is True, will prompt for username.

    Returns:
        :py:class:`networkx.Graph`: A networkx object of the desired network
        or :py:class:`ndex2.nice_cx_network.NiceCXNetwork`: If return_cx is True, returns the CX object instead of a networkx object.
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


def upload_network(G, name, template=None, networkset=None, networkset_name=None, username=None, password=None, is_cx=False):
    """ Wrapper function for uploading a network to NDEx.
    Args:
        G (networkx.Graph or ndex2.nice_cx_network.NiceCXNetwork): The network to upload.
        name (str): The name of the network.
        template (str, optional): The uuid of a template network to use as a style reference to apply to the network. Defaults to None.
        networkset (str, optional): The ID of the network set to add the network to. Defaults to None. Set to 'create' to create a new network set.
        networkset_name (str, optional): The name of the network set to create if it does not exist. Defaults to None.
        username (str, optional): Username for NDEx. If not provided, will prompt for username.
        password (str, optional): Password for NDEx. If not provided, will prompt for password.
        is_cx (bool, optional): If True, G is already a CX object. Defaults to False.
    Returns:
        str: The NDEx identifier of the uploaded network.
        or tuple: If networkset is specified, returns a tuple of the network set ID and the network ID.
    """
    if (username is None) or (password is None):
        username=getpass('Username:')
        password = getpass('Password:')
    client = ndex2.client.Ndex2(username=username, password=password)
    if not is_cx:
        G_cx = ndex2.create_nice_cx_from_networkx(G)
        G_cx.add_network_attribute('name', name)
    else:
        G_cx=G
    
    SERVER = 'http://ndexbio.org'
    if template is not None:
        # apply a style template from an existing network
        G_cx.apply_template(SERVER, template, username=username, password=password)
    id_str = G_cx.upload_to(SERVER, username=username, password=password)
    id_str = id_str.split('/')[-1]
    if networkset is None:
        return id_str
    elif networkset =='create':
        # Initiate a new network set and add the network to this set
        set_id = client.create_networkset(networkset_name, description=networkset_name)
        set_id = set_id.split('/')[-1]
        client.add_networks_to_networkset(set_id, [id_str])
        return (set_id, id_str)
    else:
        # add network to existing network set
        client.add_networks_to_networkset(networkset, [id_str])
        return id_str
    
    
def load_saved_network_nodes(indir, net_name):
    ''' Load the list of network nodes from a saved file.
    Args:
        indir (str): Directory where the network files are saved.
        net_name (str): Name of the network, used to find the corresponding nodes file.
    Returns:
        list: A list of node IDs if the file exists, otherwise None.
    '''
    if os.path.exists(os.path.join(indir, net_name+ "_nodes.txt")):
        pc_nodes = pd.read_csv(os.path.join(indir, net_name+ "_nodes.txt"), sep='\t', header=None, index_col=0).index.tolist()
        pc_nodes = [int(x) for x in pc_nodes]
        return pc_nodes
    else:
        return

def load_saved_network_degrees(indir, net_name):
    ''' Load the degree information of the network from a saved file.
    Args:
        indir (str): Directory where the network files are saved.
        net_name (str): Name of the network, used to find the corresponding degrees file.
    Returns:
        dict: A dictionary mapping node IDs to their degrees if the file exists, otherwise None.
    '''
    if os.path.exists(os.path.join(indir, net_name+ "_degrees.txt")):
        degree_map = pd.read_csv(os.path.join(indir, net_name+ "_degrees.txt"), sep='\t', header=None, index_col=0).to_dict()[1]
        return degree_map
    else:
        return

def create_saved_nodes_and_degrees(uuid, outdir, net_name, nodes=True, degrees=True, represents=False, alias='ncbigene:'):
    ''' Create and save the nodes and degrees of a network to files.
    Args:
        uuid (str): NDEx identifier of the network.
        outdir (str): Directory where the files will be saved.
        net_name (str): Name of the network, used for naming the output files.
        nodes (bool, optional): If True, saves the nodes of the network. Defaults to True.
        degrees (bool, optional): If True, saves the degrees of the network. Defaults to True.
        represents (bool, optional): If True, uses the 'represents' attribute instead of 'GeneID' for node IDs. Defaults to False.
        alias (str, optional): If provided, uses this alias to filter the 'represents' attribute. Defaults to 'ncbigene:'.
    Returns:
        None: The function saves the nodes and degrees to files in the specified directory.
    '''
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


def create_ncbi_symbol_map_from_network(uuid, outdir, net_name):
    ''' Create a mapping of Entrez IDs to gene symbols from a network. Requires the network to have 'GeneID' as a node attribute, and be indexed by symbols.
    Args:
        uuid (str): NDEx identifier of the network.
        outdir (str): Directory where the mapping file will be saved.
        net_name (str): Name of the network, used for naming the output file.
    Returns:
        None: The function saves the mapping to a file in the specified directory.
    '''
    G = load_network(uuid)
    node_map = {x[1]['GeneID']: x[0] for x in G.nodes(data=True)}
    node_df = pd.DataFrame({'Symbol':node_map}).reset_index().rename(columns={'index':'Entrez'})
    node_df.to_csv(os.path.join(outdir, f'{net_name}_node_map.txt'), sep='\t', index=False)
