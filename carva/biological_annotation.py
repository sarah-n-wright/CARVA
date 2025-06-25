import os
import pandas as pd
import argparse
from geneset_utils import *
from network_utils import *
from tqdm import tqdm

class GeneSet:
    def __init__(self, trait, indir):
        self.trait = trait
        self.cv = load_seed_genes(self.trait, 'common', indir, usecol='Gene Symbol')
        self.rv = load_seed_genes(self.trait, 'rare', indir, usecol='Gene Symbol')
        self.set_sizes = {'common':len(self.cv), 
                            'rare':len(self.rv), 
                            'overlap':len(set(self.cv).intersection(set(self.rv)))}
        



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create gene sets for a trait')
    parser.add_argument('--indir', type=str, help='Input directory')
    parser.add_argument('--outdir', type=str, help='Output directory')
    args = parser.parse_args()
    softi_dir = '/cellar/users/snwright/Git/Network_Evaluation_Tools/'
    expression_data = os.path.join(softi_dir, 'Data', 'gtex_median_processed_1.tsv.gz')
    abundance_data = ''
    citation_data = ''
    conservation_data = os.path.join(softi_dir,'StateOfTheInteractomes_Notebooks/Data/', 'gene_conservation_scores.txt')
    