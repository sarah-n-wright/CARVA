"""
Calculate gene overlap significance between rare and common variant gene sets.

This script performs hypergeometric tests to assess whether the overlap between
two gene sets (typically rare and common variants for traits) is statistically
significant given a background gene set size.

File naming convention:
    - Rare variant files: <trait>_RV.txt
    - Common variant files: <trait>_CV.txt

Required file format (two options):
    Option A (with header): Tab-delimited with 'Entrez' and 'P-value' columns
    Option B (no header): Tab-delimited with gene IDs in the first column

Usage:
    python gene_overlap.py --raretrait TRAIT1 --commontrait TRAIT2 \
        --datadir /path/to/data --outdir /path/to/output \
        --rare_th 0.05 --common_th 0.05

Output:
    Creates a TSV file: <raretrait>_<commontrait>.R_C_overlap.txt
    Columns: trait_pair, test_name, n_common, n_rare, n_intersect, background_N, p_value
"""

import argparse
import os
from scipy.stats import hypergeom
from geneset_utils import load_gene_profile


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--raretrait', type=str, help='Trait1 to evaluate')
    parser.add_argument('--commontrait', type=str, help='Trait2 to evaluate')
    parser.add_argument('--datadir', type=str, help='Path to data')
    parser.add_argument('--outdir', type=str, help='Path to output')
    parser.add_argument('--rare_th', type=float, default=1, required=False, help='P-value threshold for including a gene association in the gene set')
    parser.add_argument('--common_th', type=float, default=1, required=False, help='P-value threshold for including a gene association in the gene set')
    parser.add_argument('--min_genes', type=int, default=3, required=False, help='Minimum number of genes')
    parser.add_argument('--background_N', type=int, default=20000, required=False, help='Background gene set size')
    parser.add_argument('--test_name', type=str, help='Name to give the test')
    args = parser.parse_args()

    outfile=os.path.join(args.outdir, f'{args.raretrait}_{args.commontrait}.R_C_overlap.txt')
    rare_genes = load_gene_profile(os.path.join(args.datadir, args.raretrait + '_RV.txt'), args.rare_th, args.min_genes)
    common_genes = load_gene_profile(os.path.join(args.datadir, args.commontrait + '_CV.txt'), args.common_th, args.min_genes)
    if (rare_genes is not None) and (common_genes is not None):
        n_rare = len(rare_genes)
        n_common = len(common_genes)
        intersect = set(rare_genes).intersection(common_genes)
        n_intersect = len(intersect)
        p = hypergeom.sf(k=n_intersect-1, N=n_rare, n=n_common, M=args.background_N)
            
        with open(outfile, 'w') as f:
            f.write('\t'.join([f'{args.raretrait}_{args.commontrait}', args.test_name, str(n_common), str(n_rare), str(n_intersect), str(args.background_N), str(p)])+'\n')
        
