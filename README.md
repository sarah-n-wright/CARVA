# CARVA (*C*ommon *A*nd *R*are *V*ariant *A*nalysis)

Wright SN, Yang J, Ideker T. Common and rare genetic variants show network convergence for a majority of human traits. medRxiv (2025). DOI: [https://doi.org/10.1101/2025.06.27.25330419](https://doi.org/10.1101/2025.06.27.25330419)

## DEPENDENCIES
* python >= 3.10
* netcoloc v0.1.8 ([https://github.com/ucsd-ccbb/NetColoc/tree/quant_netcoloc](https://github.com/ucsd-ccbb/NetColoc/tree/quant_netcoloc))
* sentence-transformers 3.4.1
* ndex2 3.9.0
* neteval v0.2.2 ([https://pypi.org/project/neteval/](https://pypi.org/project/neteval/))
* networkx 2.8.8
* obonet v1.0.0
* scikit-learn v1.5.1
* R >= 4.0
* TissueEnrich ([Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/TissueEnrich.html))

## DATA GUIDE

The following data is included in the `Reference_Data/` directory:  

### GWAS Catalog Data
| File | Description |
|-----| -------|
|gwas-catalog-v1.0.3.1-studies-r2025-03-26.tsv.gz          | Raw study information for all studies in the GWAS Catalog, downloaded on March 26, 2025 |
|gwas_catalog_trait-mappings_r2025-03-26.tsv.gz            | Raw trait information for all traits in the GWAS Catalog, downloaded on March 26, 2025 |
|cleaned_gwas-catalog-v1.0.3.1-studies-r2025-03-26.tsv.gz  | Study information for filtered studies in the GWAS Catalog, with updated EFO mappings using trait information |
|gwas_catalog_Jan29_2025.txt.gz                            | Raw association data for all studies and traits, downloaded on January 29, 2025 |
|gwas_catalog_Jan29_2025.txt.cleaned.gz                    | Filtered association data |
|gwas_catalog_Jan29_2025.txt.cleaned.entrez.gz             | Filtered association data, with gene identifiers mapped to NCBI Gene IDs |

### RAVAR Data
| File | Description |
|-----| -------|
|gene_fulltable_06112024.txt.gz                      | Raw association data for all studies and traits, downloaded on June 11, 2024 |
|trait_allinfo_06112024.txt                          | Raw trait information for all traits in RAVAR, downloaded on June 11, 2024 |
|rv_study_info_cleaned_with_manual.tsv               | Manually compiled study information for all studies in RAVAR |
|rv_study_info_cleaned_with_manual_mapped_Mar28.tsv  | Study information with cleaned EFO mappings | 
|gene_fulltable_06112024.txt.entrez.gz               | Filtered association data with gene identifiers mapped to NCBI Gene IDs |

### Annotation Data
| File | Description |
|-----| -------|
|Ensembl_Feb14_2025.txt.gz               | Ensembl gene annotations, downloaded February 14, 2025 |
|IHME-GBD_2021_DATA-Prevalence.csv       | Global disease prevalence estimates, downloaded April 29, 2025 |
|gene2go.gz                              | Gene - GO mappings, downloaded November 25, 2024 |
|gnomad.v4.1.constraint_metrics.tsv.gz  | Gene constraint metrics v4.1, downloaded April 14, 2024 |
|gtex_median_processed_1.tsv.gz         | Median mRNA Expression data, generated November 25, 2024 |
|ukb31063_h2_topline.02Oct2019.tsv.gz   | SNP-heritability estimates from the UK Biobank, downlowed April 29, 2024 |

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
