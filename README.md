# CARVA (*C*ommon *A*nd *R*are *V*ariant *A*nalysis)

Wright SN, Yang J, Ideker T. Common and rare genetic variants show network convergence for a majority of human traits. bioRxiv (2025). DOI: [insert link](link)

## DEPENDENCIES
* netcoloc v0.1.8 ([https://github.com/ucsd-ccbb/NetColoc/tree/quant_netcoloc](https://github.com/ucsd-ccbb/NetColoc/tree/quant_netcoloc))
* sentence-transformers 3.4.1
* ndex2 3.9.0
* neteval v0.2.2 ([https://pypi.org/project/neteval/](https://pypi.org/project/neteval/))
* networkx 2.8.8
* obonet v1.0.0
* scikit-learn v1.5.1
* TissueEnrich ([Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/TissueEnrich.html))

## NOTEBOOK GUIDE

### `1A_DataCuration_GWAS_and_RV_processing.ipynb`
Processing of common variant associations from the GWAS Catalog and rare variant gene-level associations from RAVAR.bio. Related to Figure 1 & Supplemental Table 1. 

### `1B_DataCuration_StudySelection.ipynb`
Filtering based on study metadata and matching of common variant studies to rare variant studies for the same trait. Related to Figure 1 & Supplemental Table 1.

### `1C_DataCuration_Annotations.ipynb`
Curation of gene properties, and analysis of gene level overlap between common variant and rare variant associated genes for the same trait. Related to Figure 1. 

### `2_NetworkColocalization.ipynb`
Results of network colocalization analysis across 373 human traits. Related to Figure 2 & Supplemental Table 2. 

### `3_Reproducibility.ipynb`
Results of network colocalization analysis across 373 human traits, including results from multiple studies of the same trait, and analysis with multiple networks. Related to Figure 3, Supplemental Figure 1 & Supplemental Table 2. 

### `4A_FeaturePreparation.ipynb`
Curation of gene and trait features for each pair of common and rare variant studies. Related to Figure 4, Supplemental Figure 2, Supplemental Table 3 and Supplemental Table 4.

### `4B_FeatureAnalysis.ipynb`
Regression analysis of gene and trait features for prediction of network colocalization. Related to Figure 4, Figure 5, Supplemental Figure 2 and Supplemental Table 4. 

### `5_Examples_NeuropsychiatricTraits.ipynb`
Identification and analysis of network colocalization for neuropsychiatric traits. Related to Figure 6 & Supplemental Figure 3.

### Supplemental Notebooks
**`SX_ColocalizedNetwork_Uploads.ipynb`** - Formation of colocalized trait networks and upload to NDEx.  
**`SX_HCX_Creation.ipynb`** - Creation of hierarchical systems maps in HCX format.  
**`SX_NetColocOptimization.ipynb`** - Optimization and benchmarking of NetColoc. Related to Supplemental Figure 4 & Supplemental Table 5. 
