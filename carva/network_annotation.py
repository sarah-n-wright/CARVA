import os
import pandas as pd
import argparse
from geneset_utils import *
from network_utils import *
from tqdm import tqdm
import sys
from neteval.Timer import Timer

###################
# Gene Annotation #
###################
class GeneSet:
    def __init__(self, trait1,  indir, trait2=None, node_col='Gene Symbol'):
        if trait2 is None:
            self.trait = trait1
            self.cv = load_seed_genes(self.trait, 'common', indir, usecol=node_col)
            self.rv = load_seed_genes(self.trait, 'rare', indir, usecol=node_col)
        else:
            self.rv = load_seed_genes(trait1, 'rare', indir, usecol=node_col)
            self.cv = load_seed_genes(trait2, 'common', indir, usecol=node_col)
            self.trait = trait1 + '_' + trait2
        self.ov = set(self.cv).intersection(set(self.rv))
        self.set_sizes = {'common':len(self.cv), 
                            'rare':len(self.rv), 
                            'overlap':len(set(self.cv).intersection(set(self.rv)))}

    def get_all_stats(self, network, timer=None):
        if timer is not None:
            timer.start(f'Rare')
        rv_stats = get_subnetwork_stats(network, self, 'rare').T
        if timer is not None:
            timer.end(f'Rare')
            timer.start(f'Common')
        cv_stats = get_subnetwork_stats(network, self, 'common').T
        if timer is not None:
            timer.end(f'Common')
            timer.start(f'Overlap')
        ov_stats = get_subnetwork_stats(network, self, 'overlap').T
        if timer is not None:
            timer.end(f'Overlap')
            timer.start(f'RareCommon')
        rc_stats = get_rare_common_network_stats(network, self).T
        if timer is not None:
            timer.end(f'RareCommon')
        stats_df = pd.concat([rv_stats, cv_stats, ov_stats, rc_stats], axis=0)
        stats_df.reset_index(inplace=True, names='trait')
        return stats_df

    def get_modularity(self, network, timer=None, gamma=1):
        if timer is not None:
            timer.start('Common')
        Q_common = network.get_modularity(self.cv, gamma=gamma)
        if timer is not None:
            timer.end('Common')
            timer.start('Rare')
        Q_rare = network.get_modularity(self.rv, gamma=gamma)
        if timer is not None:
            timer.end('Rare')
            timer.start('Overlap')
        Q_overlap = network.get_modularity(self.ov, gamma=gamma)
        if timer is not None:
            timer.end('Overlap')
            timer.start('RC')
        Q_all = network.get_modularity(set(self.cv).union(self.rv), gamma=gamma)
        if timer is not None:
            timer.end('RC')
            timer.start('SubNet')
        Q_subnetwork= network.get_subnetwork_modularity(self.rv, self.cv, self.ov, gamma=gamma)
        if timer is not None:
            timer.end('SubNet')
        return {'Modularity_common':Q_common, 'Modularity_rare': Q_rare, 'Modularity_overlap':Q_overlap, 
                'Modularity_subnetwork':Q_subnetwork, 'Modularity_all': Q_all}        
    

    
class NDExNetwork:
    def __init__(self, net_input, net_name, username=None, password=None, use_lcc=True, clustering_file=None, paths_file=None, 
                 input_is_uuid=False, netcol1='Entrez_A', netcol2='Entrez_B'):
        self.net_input = net_input
        if input_is_uuid:
            self.full_G = load_network(self.net_input, username, password)
        else:
            self.full_G = nx.from_pandas_edgelist(pd.read_csv(net_input), source=netcol1, target=netcol2)
        self.net_name = net_name
        self.n_components = nx.number_connected_components(self.full_G)
        self.is_lcc = False
        if use_lcc:
            self.lcc = max(nx.connected_components(self.full_G), key=len)
            self.G = self.full_G.subgraph(self.lcc)
            self.is_lcc=True
            print(f'Using LCC with {len(self.lcc)}/{len(self.full_G.nodes())} nodes')
        else:
            self.lcc=None
            self.G = self.full_G
        self.degree_map = dict(self.G.degree())
        self.average_degree = sum(self.degree_map.values())/len(self.degree_map)
        self.density = nx.density(self.G)
        self.node_count = len(self.G.nodes)
        self.edge_count = len(self.G.edges)
        self.nodes = list(self.G.nodes)
        if (clustering_file is None) or (not os.path.exists(clustering_file)):
            self.clustering_coeffs = pd.DataFrame({'clustering_coeff': nx.clustering(self.G)})
            if clustering_file is not None:
                self.clustering_coeffs.to_csv(clustering_file)
        else:
            self.clustering_coeffs = pd.read_csv(clustering_file, index_col=0)
            if netcol1 == 'Entrez_A':
                self.clustering_coeffs.index = self.clustering_coeffs.index.astype(int)
            self.clustering_coeffs.columns = ['clustering_coeff']
        if (paths_file is not None) and (os.path.exists(paths_file)):
            self.path_lengths = pd.read_csv(paths_file, index_col=0)
            self.average_path_lengths = self.path_lengths.sum().sum()/((len(self.path_lengths)**2 - len(self.path_lengths))) # excludes self paths (which are 0 by definition)
            if netcol1 == 'Entrez_A':
                self.path_lengths.index = self.path_lengths.index.astype(int)
                self.path_lengths.columns = [int(x) for x in self.path_lengths.columns]
        else:
            self.path_lengths = None
            self.average_path_lengths= np.nan
        self.clustering_coefficient = self.clustering_coeffs.clustering_coeff.mean()
            

    def get_average_clustering(self, genelist):
        present_nodes = [g for g in genelist if g in self.clustering_coeffs.index.values]
        try:
            return self.clustering_coeffs.loc[present_nodes].clustering_coeff.mean()
        except:
            print(present_nodes)
            print(self.clustering_coeffs.index.values[0:5])
            print(isinstance(present_nodes[0], int))
            print(isinstance(self.clustering_coeffs.index.values[0], np.int64))
            return self.clustering_coeffs.loc[present_nodes].clustering_coeff.mean()
    
    def get_average_shortest_path(self, genelist):
        if self.path_lengths is not None:
            try:
                present_nodes = [g for g in genelist if g in self.path_lengths.index.values]
                return self.path_lengths.loc[present_nodes, present_nodes].sum().sum()/(len(present_nodes)**2 - len(present_nodes))
            except:
                present_nodes = [g for g in genelist if g in self.path_lengths.index.values]
                print(present_nodes)
                print(isinstance(present_nodes[0], int))
                print(isinstance(self.path_lengths.index.values[0], np.int64))
                print(self.path_lengths.index.values[0:5])
                return self.path_lengths.loc[present_nodes, present_nodes].sum().sum()/(len(present_nodes)**2 - len(present_nodes))
        else:
            return np.nan
    
    
    def get_nodes_in_network(self, genelist):
        present_nodes = [gene for gene in genelist if gene in self.nodes]
        if len(present_nodes) < len(genelist):
            print(f'{len(present_nodes)}/{len(genelist)} genes found in network')
        return present_nodes
    
    def get_node_degrees(self, genelist):
        return {gene: self.degree_map[gene] for gene in genelist if gene in self.degree_map}
    
    
    def export_global_network_stats(self):
        stats = {'uuid': self.uuid,
                    'lcc': self.is_lcc,
                    'node_count': self.node_count,
                    'edge_count': self.edge_count,
                    'density': self.density,
                    'clustering_coefficient': self.clustering_coefficient,
                    'average_degree': self.average_degree}
        return pd.DataFrame({self.net_name: stats})
    
    def get_subnetwork(self, genelist):
        subgraph = self.G.subgraph(genelist)
        return subgraph

    def get_modularity(self, genelist, gamma=1):
        present_nodes = self.get_nodes_in_network(genelist)
        communities = [present_nodes, [node for node in self.G if node not in present_nodes]]
        return nx.community.modularity(self.G, communities, resolution=gamma)

    def get_subnetwork_modularity(self, rare_list, common_list, overlap_list, gamma=1):
        partition_o = self.get_nodes_in_network(overlap_list)
        partition_r = [x for x in self.get_nodes_in_network(rare_list) if x not in partition_o]
        partition_c = [x for x in self.get_nodes_in_network(common_list) if x not in partition_o]
        communities = []
        for partition in [partition_r, partition_c, partition_o]:
            if len(partition) > 0:
                communities.append(partition)
        all_nodes = set(partition_r + partition_c + partition_o)
        subG = self.get_subnetwork(all_nodes)
        if len(subG.edges) > 0:
            return nx.community.modularity(subG, communities, resolution=gamma)
        else:
            print(f'No edges in subnetwork. Returning modularity=0')
            return 0


    def get_assortativity(self, genelist):
        mx_graph = self.G.copy()
        attributes = {node: 1 if node in genelist else 0 for node in mx_graph.nodes}
        nx.set_node_attributes(mx_graph, attributes, 'InGeneSet')
        return nx.attribute_assortativity_coefficient(mx_graph, 'InGeneSet')
    
    def get_subgraph_clustering(self, genelist):
        return nx.average_clustering(self.G.subgraph(genelist))

def get_subnetwork_stats(network, geneset, rare_or_common):
    if rare_or_common == 'rare':
        genes_rv = network.get_nodes_in_network(geneset.rv)
        if len(genes_rv) == 0:
            return pd.DataFrame({geneset.trait: {'varset': 'rare',
                                            'subnetwork_density': 0,
                                            'n_components': 0,
                                            'average_degree': 0,
                                            'subnetwork_average_degree': 0,
                                            'clustering':0, 
                                            'subnetwork_clustering': 0,
                                            'assortativity': 0,
                                            'average_path': 0}})
        sub_G = network.get_subnetwork(genes_rv)

    elif rare_or_common == 'common':
        genes_cv = network.get_nodes_in_network(geneset.cv)
        if len(genes_cv) == 0:
            return pd.DataFrame({geneset.trait: {'varset': 'common',
                                            'subnetwork_density': 0,
                                            'n_components': 0,
                                            'average_degree': 0,
                                            'subnetwork_average_degree': 0,
                                            'clustering':0,
                                            'subnetwork_clustering': 0,
                                            'assortativity': 0,
                                            'average_path': 0}})
        sub_G = network.get_subnetwork(genes_cv)
    elif rare_or_common == 'overlap':
        genes_ov = list(set(network.get_nodes_in_network(geneset.cv)).intersection(set(network.get_nodes_in_network(geneset.rv))))
        if len(genes_ov) == 0:
            return pd.DataFrame({geneset.trait: {'varset': 'overlap',
                                            'subnetwork_density': 0,
                                            'n_components': 0,
                                            'average_degree': 0,
                                            'subnetwork_average_degree': 0,
                                            'clustering':0,
                                            'subnetwork_clustering': 0,
                                            'assortativity': 0,
                                            'average_path': 0}})
        sub_G = network.get_subnetwork(genes_ov)
    sub_density = nx.density(sub_G)
    sub_n_components = nx.number_connected_components(sub_G)
    degrees = network.get_node_degrees(sub_G.nodes)
    average_degree = sum(degrees.values())/len(sub_G.nodes)
    sub_average_degree = sum(dict(sub_G.degree()).values())/len(sub_G.nodes)
    clustering_coefficient = nx.average_clustering(sub_G)
    stats = {'varset': rare_or_common, 
                'subnetwork_density': sub_density,
                'n_components': sub_n_components,
                'average_degree': average_degree,
                'subnetwork_average_degree': sub_average_degree,
                'clustering': network.get_average_clustering(sub_G.nodes),
                'subnetwork_clustering': clustering_coefficient,
                'assortativity': network.get_assortativity(sub_G.nodes),
                'average_path': network.get_average_shortest_path(sub_G.nodes)}
    return pd.DataFrame({geneset.trait: stats})

def get_rare_common_network_stats(network, geneset):
    r_genes = set(network.get_nodes_in_network(geneset.rv))
    c_genes = set(network.get_nodes_in_network(geneset.cv))
    if len(r_genes) == 0 or len(c_genes) == 0:
        return pd.DataFrame({geneset.trait: { 'varset': 'rare_common',
                                                'subnetwork_density': 0,
                                                'n_components': 0,
                                                'average_degree': 0,
                                                'subnetwork_average_degree': 0,
                                                'clustering': 0,
                                                'subnetwork_clustering': 0,
                                                'subnetwork_assortativity': 0,
                                                'average_path': 0}})
    o_genes = r_genes.intersection(c_genes)
    all_genes = r_genes.union(c_genes)
    sub_G = network.get_subnetwork(list(all_genes))
    sub_density = nx.density(sub_G)
    sub_n_components = nx.number_connected_components(sub_G)
    sub_average_degree = sum(dict(sub_G.degree()).values())/len(sub_G.nodes)
    average_degree = sum(network.get_node_degrees(all_genes).values())/len(all_genes)
    clustering_coefficient = nx.average_clustering(sub_G)
    sub_assortativity = rare_common_assortativity(sub_G, r_genes, c_genes, o_genes)
    stats = {'varset': 'rare_common',
                'subnetwork_density': sub_density,
                'n_components': sub_n_components,
                'average_degree': average_degree,
                'subnetwork_average_degree': sub_average_degree,
                'clustering': network.get_average_clustering(sub_G.nodes),
                'subnetwork_clustering': clustering_coefficient,
                'assortativity': network.get_assortativity(sub_G.nodes), 
                'subnetwork_assortativity': sub_assortativity,
                'average_path': network.get_average_shortest_path(sub_G.nodes)}
    return pd.DataFrame({geneset.trait: stats})

def rare_common_assortativity(network, r_genes, c_genes, o_genes):
    mx_graph = network.copy()
    if len(o_genes) == 0:
        attributes = {node: 1 if node in r_genes else 0 for node in mx_graph.nodes}
        nx.set_node_attributes(mx_graph, attributes, 'InGeneSet')
        return nx.attribute_assortativity_coefficient(mx_graph, 'InGeneSet')
    else:
        # Assess Rare+overlap vs Common-only , and Common+overlap vs Rare only
        att_RO_C = {node: 1 if (node in o_genes) or (node in r_genes) else 0 for node in mx_graph.nodes}
        att_CO_R = {node: 1 if (node in o_genes) or (node in c_genes) else 0 for node in mx_graph.nodes}
        nx.set_node_attributes(mx_graph, att_RO_C, 'InGeneSet_RO_C')
        nx.set_node_attributes(mx_graph, att_CO_R, 'InGeneSet_CO_R')
        assortativity_RO_C = nx.attribute_assortativity_coefficient(mx_graph, 'InGeneSet_RO_C')
        assortativity_CO_R = nx.attribute_assortativity_coefficient(mx_graph, 'InGeneSet_CO_R')
        return (assortativity_RO_C, assortativity_CO_R)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, default = '',
                        help='Path to inputs')
    parser.add_argument('--outdir', type=str, default = '', help='Path to outputs')
    parser.add_argument('--uuid', type=str, default=None,
                        help='UUID of network to analyze')
    parser.add_argument('--netfile', type=str, default=None)
    parser.add_argument('--geneset_list_file', type=str, default=None)
    parser.add_argument('--geneset_list_file2', type=str, default=None)
    parser.add_argument('--clustering_file', type=str, default=None)
    parser.add_argument('--paths_file', type=str, default=None)
    parser.add_argument('--net_name', type=str, default='default')
    parser.add_argument('--update', action='store_true', help='Force update of already calculated values.')
    args = parser.parse_args()
    t = Timer()
    t.start('Load Network')
    
    assert (args.uuid is not None) or (args.netfile is not None), 'Either UUID or network file must be specified'
    
    if args.netfile is not None:
        G_nd = NDExNetwork(args.netfile, args.net_name, clustering_file=args.clustering_file, paths_file=args.paths_file, input_is_uuid=False, netcol1='Entrez_A', netcol2='Entrez_B')
    else:
        try:
            G_nd = NDExNetwork(args.uuid, args.net_name, clustering_file=args.clustering_file, paths_file=args.paths_file, input_is_uuid=True)
        except:
            try:
                G_nd = NDExNetwork(args.uuid, args.net_name, clustering_file=args.clustering_file, paths_file=args.paths_file, input_is_uuid=True)
            except:
                G_nd = NDExNetwork(args.uuid, args.net_name, clustering_file=args.clustering_file, paths_file=args.paths_file, input_is_uuid=True)
    t.end('Load Network')
    
    if args.geneset_list_file2 is None:
        with open(os.path.join(args.indir, args.geneset_list_file), 'r') as f:
            all_genesets = f.read().splitlines()
        all_genesets2 = [None for _ in all_genesets]
        outprefs = [g for g in all_genesets]
    else:
        with open(os.path.join(args.indir, args.geneset_list_file), 'r') as f:
            all_genesets = f.read().splitlines()
        with open(os.path.join(args.indir, args.geneset_list_file2), 'r') as f:
            all_genesets2 = f.read().splitlines()
            
        outprefs = [all_genesets[i] + '_' + all_genesets2[i] for i in range(len(all_genesets))]
    print('OUTPREFs', outprefs)
    results = []
    
    if True:
        for i, geneset in tqdm(enumerate(all_genesets)):
            if (not args.update) and (os.path.exists(os.path.join(args.outdir, f'network_stats_{args.net_name}_{args.geneset_list_file}.{outprefs[i]}'))):
                # check if already calculated, skip if true
                continue
            t.start(f'GeneSet {outprefs[i]}') 
            gs = GeneSet(geneset, args.indir, trait2=all_genesets2[i], node_col='Entrez')
            geneset_results = gs.get_all_stats(G_nd, timer=t)
            results.append(geneset_results)
            geneset_results.to_csv(os.path.join(args.outdir, f'network_stats_{args.net_name}_{args.geneset_list_file}.{outprefs[i]}'), sep='\t', index=False)
            t.end(f'GeneSet {outprefs[i]}')
        
        all_results = pd.concat(results, axis=0)
        all_results['Network'] = args.net_name
        all_results.to_csv(os.path.join(args.outdir, f'network_stats_{args.net_name}_{args.geneset_list_file}'), sep='\t', index=False)
        t.print_all_times()
    if True:
        results = {}
        for i, geneset in tqdm(enumerate(all_genesets)):
            if (not args.update) and (os.path.join(args.outdir, f'network_modularities_{args.net_name}_{args.geneset_list_file}')):
                # check if already calculated, skip if true
                continue
            t.start(f'GeneSet {outprefs[i]}') 
            gs = GeneSet(geneset, args.indir, trait2=all_genesets2[i], node_col = 'Entrez')
            x = gs.get_modularity(G_nd, timer=t, gamma=1)
            results[outprefs[i]] = x
        results_df = pd.DataFrame(results).T
        results_df.to_csv(os.path.join(args.outdir, f'network_modularities_{args.net_name}_{args.geneset_list_file}'), sep='\t', index=False)