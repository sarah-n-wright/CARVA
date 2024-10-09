import pandas as pd
import networkx as nx
import os
import ndex2

datadir = '/cellar/users/snwright/Data/NetColocTest/'

ndex_server='public.ndexbio.org'

ndex_user=None
ndex_password=None
G_overlap_cx = ndex2.create_nice_cx_from_server(
        ndex_server, 
        username=ndex_user, 
        password=ndex_password, 
        uuid='d73d6357-e87b-11ee-9621-005056ae23aa')
G_overlap = G_overlap_cx.to_networkx()
print('number of nodes:')
print(len(G_overlap.nodes))
print('\nnumber of edges:')
print(len(G_overlap.edges))

path_lengths = nx.all_pairs_shortest_path_length(G_overlap)
pd.DataFrame.from_dict(dict(path_lengths)).to_csv(os.path.join(datadir, 'Reference', 'path_lengths1.csv'))
