# Notebook Manifest

---
## Input Data

#### `RAVAR_GWASCat.ipynb`

* Load and process EFO ontology
  * *In*: `'http://www.ebi.ac.uk/efo/efo.obo'` 
* Load and process RAVAR data
  * *In*: `RAVAR/gene_fulltable_06112024.txt`
  * *Out*: 
* Load and process GWAS catalog data
  * *In*: `GWASCatalog/gwas_cat_download_Jan22.txt`
  * *Out*: `GWASCatalog/gwas_cat_download_Jan22.txt.cleanedJun17`
* Convert all genes to Entrez
  * *In*: `ravar['Ensembl ID']`, `gwas_genes['SNP_GENE_ID]`
  * *Out*: `<file>.entrez`
* Identify genes in both input data sets
* Descriptive statistics of input data
    * Distributions of P-values
    * Number of genes per trait
* Export gene profiles
  * *In*: dta with converted gene indentifiers
  * *Out*: `RareCommon/inputs/<trait>_CV.txt`, `RareCommon/inputs/<trait>_RV.txt`
* Trait overlap
  * *Out*: `RareCommon/inputs/overlap_traits_June27.txt'
  * Cumulative number of traits with > n genes in R, C and RC
* Trait Information
    * Classification of traits into parent categories

**Figures:**
 * Overall venn diagrams
 * Distributions of trait genes
 * Trait Classification bar charts
---
#### `go_genesets.ipynb`

Note: Utilizes `neteval` code for managing GO Data and was reused for EGAD go genesets

* Load GO data
  * *In*: `https://current.geneontology.org/ontology/go.obo`, NCBI associationa 
* Subset to terms within each branch with between n and N genes
    * *Out*: `NetColocTest/Reference/go.genesets`
* Ontology distance between terms
    * Indentify the minimum distance between terms as a measure of how similar they are
    * *Out*: `NetColocTest/inputs/trait_pairs_cross_GO.txt` - Pairs of traits that are distal within the onotology. 
* Generation of GO gene sets
    * Assignment to R/C is arbitrary. 
    * *Out*: `NetColocTest/inputs/GO/<trait>_[RV|CV].txt`
    * #TODO figure out exactly what the hypothesis here was

---
#### `debug_pcnet2.ipynb`

Note exactly sure what the point with this was. 
* Generates and debugs degree distribution of PCNet2
    * *In*: `NetColocTest/Reference/pcnet_2_0_degrees.txt`
    * *Out*: `RareCommon/inputs/pcnet2_0_degree.txt`
* Assign degrees to bins
* Loads pilot results, and identifies within trait pairs
    * *In*: `RareCommon/outputs/netcoloc_pilot_results.tsv`
    * *Out*: `RareCommon/outputs/pilot_netcoloc_results_within.tsv`

---
## Methods Development

#### #TBD

---
## Results & Figures

#### `Gene_Overlap.ipynb`

*Env*: CARVA

Analysis gene set overlaps
* Load high throughput results of statistical overlap between R & C
    * *In*: `RareCommon/outputs/overlap_results.<group>.txt
    * Variation in the background used for statistical testing, as well as whether genes were first subset to significant associations
* Descriptive statistics of overlap results
* Justification for looking beyond the overlap
* Look at overlap with different identifiers

**Figures:**
* Number of genes per gene set significant vs all
* Number of overlapping genes
* Scatter plots of the number of genes in each group as a function of parameters
    * bar plots of the same
* Box plots of the p value of overlap
    * Most significant overlaps with sig_only
* **Stacked bar plots for trait associations**
* **Venn diagrams of differing overlaps**
----
#### `GenesetFeatures.ipynb`

*Env*: CARVA

* Load the genesets (go.genesets)
* Testing of all shortest paths
* Testing of loading PCNet2
* Average degree of genesets

---

#### `NetColocTesting.ipynb`

*Env*: CARVA

Assessment of results from simulation studies with GO genesets

* Functions for loading and visualizing results
* **Remove** results
    * *In*: `load_netcoloc_results('overlap_only_remove.txt')`
* **Bin** results
    * *In*: `load_netcoloc_results('overlap_only_bin.txt')`
* Constant set sizes
    * *In*: `load_netcoloc_results('overlap_constant_[remove|bin].txt')`
* Linear regressions of overall trends
* Repeats test
    * *In*: `load_netcoloc_results('repeat_test.txt')`
* GO Relevance (dilution test)
    * *In*: `load_netcoloc_results('all_relevance_only.txt')`
* Cross-term results
    * *In*: `NetColocTest/outputs/cross_go_results.txt` , `NetColocTest/outputs/gene_overlap_cross_trait.txt`

**Figures:**

Many exploratory figures of the relationships between results and gene set overlap etc.

---

#### `Netcoloc.ipynb`

Initial assessment of the spectrum of results

* Load initial pilot results
    * *In*: `pilot_netcoloc_within_v0.1.txt`
    * Overall distributions
* Combine with overlap results
    * *In*: `overlap_results.defaults_9k.txt`
* Correlation between NetColoc and Overlap
* Identification of most extreme traits

**Figures:**
* **Pilot spectrum plot**
* Exploratory distributions
* Comparisons of overlap and NPS

---