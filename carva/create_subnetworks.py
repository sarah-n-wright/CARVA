import pandas as pd
import networkx as nx
import argparse
import os
from ndex2.cx2 import RawCX2NetworkFactory
from network_utils import *
from geneset_utils import *

def load_z(trait, datadir, transform, norm, rorc):
    """ Load z-scores for a given trait and return as a DataFrame.
    Args:
        trait (str): Trait name to load z-scores for.
        datadir (str): Directory containing the z-score files.
        transform (str): Transformation applied to the z-scores.
        norm (str): Normalization applied to the z-scores.
        rorc (str): 'R' for rare or 'C' for common.
    Returns:
        pd.DataFrame: DataFrame containing z-scores indexed by Entrez IDs.
    """
    assert rorc in ['R', 'C'], "rorc must be either 'R' or 'C'"
    assert transform in ['neglog10', 'log10'], "transform must be either 'neglog10' or 'log10'"
    assert norm in ['sum', 'max'], "norm must be either 'sum' or 'max'"
    df = pd.read_csv(os.path.join(datadir, f'{trait}_z_{rorc}V_q_{transform}_{norm}.tsv'), sep='\t', header=None,
                    index_col=0, names=['z']).dropna()
    return df

def load_genelists(trait, genelist_dir, rorc):
    """ Load gene lists for a given trait and return as a DataFrame.
    Args:
        trait (str): Trait name to load gene lists for.
        genelist_dir (str): Directory containing the gene list files.
        rorc (str): 'R' for rare or 'C' for common.
    Returns:
        pd.DataFrame: DataFrame containing gene lists with Entrez IDs, p-values, and
        gene symbols.
    """
    assert rorc in ['R', 'C'], "rorc must be either 'R' or 'C'"
    # Load the gene list file and rename columns
    df = pd.read_csv(os.path.join(genelist_dir, f'{trait}_{rorc}V.txt'), sep='\t')
    df = df.rename(columns={'P-value': f'pval_{rorc}', 'Gene Symbol': f'symbol_{rorc}'})
    df = df.sort_values(by=f'pval_{rorc}', ascending=True)
    df  =df.drop_duplicates(subset=['Entrez'])
    return df.loc[:, ('Entrez', f'pval_{rorc}', f'symbol_{rorc}')]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--network_uuid', type=str, required=True, help='UUID of the network to use')
    parser.add_argument('--trait_list_file', type=str, required=True, help='File containing the list of traits to process, with columns TraitR and TraitC')
    parser.add_argument('--z_dir', type=str, required=True, help='Directory containing the z-score files for traits')
    parser.add_argument('--genelist_dir', type=str, required=True, help='Directory containing the gene lists for traits')
    parser.add_argument('--outputdir', type=str, required=True, help='Directory to save the output subnetworks')
    parser.add_argument('--zth', default=1, type=float, help='Individual Z-score threshold for coloc genes')
    parser.add_argument('--zzth', default=3, type=float, help='Combined Z-score threshold for coloc genes')
    parser.add_argument('--use_cx2', action='store_true', help='Use CX2 format for output. Otherwise output will be in TSV format.')
    args = parser.parse_args()

    # Load the underlying network
    G = load_network(args.network_uuid)
    mapping = {x[0]: int(x[1]['GeneID']) for x in G.nodes(data=True)}
    H = nx.relabel_nodes(G, mapping, copy=True)
    # make sure the output directory exists
    os.makedirs(args.outputdir, exist_ok=True)

    # Load the trait list and iterate over all traits
    trait_info = pd.read_csv(args.trait_list_file, sep='\t')
    trait_list = [x for x in zip(trait_info.TraitR.values, trait_info.TraitC.values)]
    for trait in trait_list:
        # load the z scores and identify the colocalized genes
        z_rare = load_z(trait[0], args.z_dir, transform='neglog10', norm='sum', rorc='R').dropna()
        z_common = load_z(trait[1], args.z_dir, transform='neglog10', norm='sum', rorc='C').dropna()
        z_df = z_rare.join(z_common, how='inner', lsuffix='_R', rsuffix='_C')
        z_df['Z_coloc'] = z_df['z_R'] * z_df['z_C']
        z_df = z_df.reset_index()
        z_df = z_df.rename(columns={'index':'Entrez'})
        z_df['Entrez'] = z_df['Entrez'].astype(int)
        if any(z_df['Entrez'].value_counts() > 1 ): # flag any duplicated Entrez IDS
            print('DUPLICATED ENTREZ 1!!!')
            counts = z_df.Entrez.value_counts()
            dupes = counts[counts>1].index.values
            print(z_rare.loc[dupes])
            print(z_common.loc[dupes])
            print(z_df.loc[z_df.Entrez.isin(dupes)])

        # load the lists of input genes 
        r_genes = load_genelists(trait[0], args.genelist_dir, 'R')
        c_genes = load_genelists(trait[1], args.genelist_dir, 'C')
        # Define node attributes
        all_genes = r_genes.merge(c_genes, on=['Entrez'], how='outer')
        all_genes['Entrez'] = all_genes.Entrez.astype(int)
        all_genes['input_gene'] = ~all_genes.pval_R.isna() | ~all_genes.pval_C.isna() 
        all_genes['rare'] = ~all_genes.pval_R.isna()
        all_genes['common'] = ~all_genes.pval_C.isna()
        all_genes['shared'] = all_genes['common'] & all_genes['rare']
        all_genes['gene_class'] = all_genes.apply(lambda x: 'shared' if (x['common'] & x['rare']) else 'common' if x['common'] else 'rare', axis=1)
        z_df = z_df.merge(all_genes, on=['Entrez'], how='left')
        z_df['input_gene'].fillna(False, inplace=True)
        z_df['gene_class'].fillna('Network', inplace=True)
        z_df = z_df.assign(HGNC=z_df.Entrez.map(map_genes_using_network_map(z_df.Entrez.values)))
        if any(z_df['Entrez'].value_counts() > 1 ):
            print('DUPLICATED ENTREZ 2!!!')
            counts = z_df.Entrez.value_counts()
            dupes = counts[counts>1].index.values
            print(z_rare.loc[dupes])
            print(z_common.loc[dupes])
            print(all_genes.loc[all_genes.Entrez.isin(dupes)])
            print(z_df.loc[z_df.Entrez.isin(dupes)])
        # Identify all colocalized genes
        z_df['coloc_gene'] = z_df.apply(lambda x: 1 if x['Z_coloc'] > args.zzth and x['z_R'] > args.zth and x['z_C'] > args.zth else 0, axis=1)
        all_df = z_df[(z_df['coloc_gene']==1) | (z_df['input_gene'])]

        # Create a subgraph of the network
        G_sub_all = nx.subgraph(H, all_df['Entrez'].tolist())
        if len([x for x in G_sub_all.nodes()]) == 0:
            print([x for x in H.nodes()][0:10])
            print(all_df['Entrez'].tolist()[0:10])
        G_sub_inputs = nx.subgraph(H, z_df[z_df['input_gene']]['Entrez'].tolist())
        # add node_attributes to the subgraph
        for Gout in [G_sub_all, G_sub_inputs]:
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['symbol_R'].to_dict(), name='symbol_R')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['symbol_C'].to_dict(), name='symbol_C')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['pval_R'].to_dict(), name='pval_R')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['pval_C'].to_dict(), name='pval_C')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['z_R'].to_dict(), name='z_R')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['z_C'].to_dict(), name='z_C')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['Z_coloc'].to_dict(), name='Z_coloc')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['gene_class'].to_dict(), name='gene_class')
            nx.set_node_attributes(Gout, all_df.set_index('Entrez')['coloc_gene'].to_dict(), name='coloc_gene')
        if args.use_cx2:
            # export to cx2
            factory = RawCX2NetworkFactory()
            cx2_all = factory.get_cx2network(G_sub_all)
            cx2_inputs = factory.get_cx2network(G_sub_inputs)
        else:
            # save the edgelists as pandas dataframes
            G_sub_all_df = nx.to_pandas_edgelist(G_sub_all)
            G_sub_inputs_df = nx.to_pandas_edgelist(G_sub_inputs)
            # save the edgelists as tsv files
            G_sub_all_df.to_csv(os.path.join(args.outputdir, f'{"_".join(trait)}_subnetwork_all.tsv'), sep='\t', index=False)
            G_sub_inputs_df.to_csv(os.path.join(args.outputdir, f'{"_".join(trait)}_subnetwork_inputs.tsv'), sep='\t', index=False)
            # save the node attributes as a tsv file
            all_df.to_csv(os.path.join(args.outputdir, f'{"_".join(trait)}_subnetwork_all_node_attributes.tsv'), sep='\t', index=False)
            z_df[z_df['input_gene']].to_csv(os.path.join(args.outputdir, f'{"_".join(trait)}_subnetwork_inputs_node_attributes.tsv'), sep='\t', index=False)
        print(f'{"_".join(trait)} done')














