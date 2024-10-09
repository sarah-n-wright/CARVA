import pandas as pd
import numpy as np
from scipy.stats import hypergeom
import argparse

#TODO deal with the duplicate definition of datadir between here and the sh script

def load_common_seed_genes(filepath, trait='',outdir="", pval=0.05):
    """
    :param filepath:
    :param interactome_nodes:
    :param trait:
    """
    all_scores = pd.read_csv(filepath, sep="\s+")
    # subset to genes in the interactome
    #all_scores = all_scores.loc[list(np.intersect1d(all_scores.index.tolist(), interactome_nodes))]
    # Calculate bonferroni corrected pvalue (alpha=0.05)
    # Get significant genes
    print(all_scores.sort_values(by='pvalue').head())
    seeds = all_scores[all_scores['pvalue'] < pval]
    seeds = seeds.sort_values(by="pvalue")
    seeds.loc[:, ("gene_symbol")].to_csv(outdir+"common_seeds_"+str(trait)+".txt", sep="\t", header=False, index=False)
    return seeds.gene_symbol.tolist()


def calculate_overlap(seeds1, seeds2, M=18820, trait="", max_seeds=18820):
    if len(seeds1) > max_seeds:
        seeds1 = seeds1[:500]
    if len(seeds2) > max_seeds:
        seeds2 = seeds2[:500]
    seeds1 = set(seeds1)
    seeds2 = set(seeds2)
    hyper = hypergeom(M=M, n=len(seeds1), N=len(seeds2))
    intersect = seeds1.intersection(seeds2)
    p_intersect = hyper.sf(len(intersect))
    return intersect, p_intersect


def get_rare_seeds(rareData, trait, pval_cutoff=1, outdir="", binary=0):
    if binary == 1:
        binary_map = pd.read_csv("binary_pheno_map.txt", header=None, names=["Field", "FullName"])
    
    data = pd.read_csv(rareData, sep="\t")
    if binary == 1:
        data.columns = ["Gene", "Pheno", "p-value"]
        data = data[data.Pheno.isin(binary_map.FullName)]
        data = data.merge(binary_map, left_on="Pheno", right_on="FullName")
    if type(trait) != type(data.Field.dtype):
        if data.Field.dtype == int:
            trait = int(trait)
        elif data.Field.dtype == object:
            if type(trait) != str:
                trait = str(trait)
                
    data = data.loc[data["Field"] == trait]
    data = data.loc[data["p-value"] < pval_cutoff]
    data = data.drop_duplicates(subset="Gene")
    data.Gene.to_csv(outdir+"rare_seeds_"+str(trait)+".txt", sep="\t", index=False, header=False)
    seeds = data.Gene.tolist()
    return seeds


if __name__=="__main__":
    datadir="/cellar/users/snwright/Data/RareCommon/"
    parser = argparse.ArgumentParser(description="Direct_overlap")
    parser.add_argument('common_file', metavar='c', type=str, help="File path to output of pascal for common variants")
    parser.add_argument('rare_file', metavar="r", type=str, help="File path for rare variant genes")
    parser.add_argument('traitc', metavar="tc",type=str, help="name of trait")
    parser.add_argument('traitr', metavar="tr",type=str, help="name of trait")
    parser.add_argument('binary', metavar="b", type=int, help="binary trait")
    parser.add_argument('out_file', metavar="o", type=str, help="File for appending result")
    args = parser.parse_args()
    if args.binary == 1:
        common_th=0.000001
    else:
        common_th=0.00000001
    common_seeds = load_common_seed_genes(args.common_file, trait=args.traitc, outdir=datadir, pval=common_th)
    rare_seeds = get_rare_seeds(args.rare_file, trait=args.traitr, outdir=datadir, pval_cutoff=0.00001, binary=args.binary)
    overlap, p_overlap = calculate_overlap(common_seeds, rare_seeds)
    with open(args.out_file, 'a') as f:
        f.write("\t".join([args.traitc, args.traitr, str(len(common_seeds)), str(len(rare_seeds)),
        str(len(overlap)), str(p_overlap)])+"\n")