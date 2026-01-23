# CARVA (*C*ommon *A*nd *R*are *V*ariant *A*nalysis)

### OVERVIEW

While both common and rare variants contribute to the genetic etiology of complex traits, whether their impacts manifest through the same effector genes and molecular mechanisms is not well understood. With the code in this repository we perform a systematic evaluation of the network convergence of common and rare variants, and assess the factors that drive this convergence.

### CITATION

Wright SN, Yang J, Ideker T. Common and rare genetic variants show network convergence for a majority of human traits. medRxiv (2025). DOI: [https://doi.org/10.1101/2025.06.27.25330419](https://doi.org/10.1101/2025.06.27.25330419)

## INSTALLATION

1. Clone this repository:
```bash
git clone https://github.com/sarah-n-wright/CARVA.git
cd CARVA
```

2. Create a conda environment and install dependencies:
```bash
conda create -n carva python=3.10
conda activate carva
pip install netcoloc==1.1.0a1 neteval==0.2.3a1 sentence-transformers==3.4.1 \
    scikit-learn==1.5.1
conda install -c conda-forge r-base==4.2
```

Install R dependencies within R (necessary for performing TSEA):
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TissueEnrich")
```

### DEPENDENCIES
* python >= 3.10
* netcoloc v1.1.0a1 ([https://pypi.org/project/netcoloc/1.1.0a1/](https://pypi.org/project/netcoloc/1.1.0a1/))
* sentence-transformers 3.4.1
* neteval v0.2.3a1 ([https://pypi.org/project/neteval/0.2.3a1/](https://pypi.org/project/neteval/0.2.3a1/))
* scikit-learn v1.5.1
* R >= 4.0
* TissueEnrich ([Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/TissueEnrich.html))


## QUICK START

* `CARVA/Notebooks/` contains all Jupyter notebooks necessary for reproducing manuscript analysis and figures. Note that some components use pre-computed results to allow these notebooks to be run on a standard computer. See **Notebook Descriptions** for further details.

* `CARVA/carva/` contains Python scripts for carrying out the common and rare variant NetColoc analysis, and modules for data processing and analysis.

* `CARVA/run/` contains scripts for implementing the full common and rare variant NetColoc analysis. Running these scripts with a large network (such as PCNet 2.0) requires access to a high-performance computing environment. See `CARVA/run/README.md` for detailed usage.

## DATA GUIDE

Networks are available at [NDEx](https://ndexbio.org):

* PCNet 2.0: [https://doi.org/10.18119/N9JP5J](https://doi.org/10.18119/N9JP5J)
* All trait-specific colocalized networks [Network Set](https://www.ndexbio.org/#/networkset/287cafe2-1645-11f0-9806-005056ae3c32?accesskey=66e87e068c4c4709d4a89604558057f386142cc100b4eed1802a04474d2c7fa4)

Expanded View Datasets EV1-6 are available in the `outputs/` directory:

* **DatasetEV1** - All input common and rare studies, properties, and gene sets.
* **DatasetEV2** - All overlap and network colocalization results.
* **DatasetEV3** - SNP heritability and disease prevalence estimates.
* **DatasetEV4** - Trait and gene set features for all common-rare study combinations.
* **DatasetEV5** - NetColoc optimization results.
* **DatasetEV6** - NetColoc benchmarking results.

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

### Annotation & Network Data
| File | Description |
|-----| -------|
|Ensembl_Feb14_2025.txt.gz               | Ensembl gene annotations, downloaded February 14, 2025 |
|IHME-GBD_2021_DATA-Prevalence.csv       | Global disease prevalence estimates, downloaded April 29, 2025 |
|gene2go.gz                              | Gene - GO mappings, downloaded November 25, 2024 |
|gnomad.v4.1.constraint_metrics.tsv.gz  | Gene constraint metrics v4.1, downloaded April 14, 2024 |
|gtex_median_processed_1.tsv.gz         | Median mRNA Expression data, generated November 25, 2024 |
|ukb31063_h2_topline.02Oct2019.tsv.gz   | SNP-heritability estimates from the UK Biobank, downloaded April 29, 2024 |
|s_het_estimates.genebayes.tsv          | Selective constraint estimates from GeneBayes |
|pcnet2_0_nodes.txt                     | Node list for PCNet 2.0 network used for analysis |
|pcnet2_0_node_map.txt                  | Gene identifier mappings for PCNet 2.0 network nodes |

## NOTEBOOK GUIDE

### `1A_DataCuration_GWAS_and_RV_processing.ipynb`
Processing of common variant associations from the GWAS Catalog and rare variant gene-level associations from RAVAR.bio. Related to Figure 1 and Dataset EV1. 

### `1B_DataCuration_StudySelection.ipynb`
Filtering based on study metadata and matching of common variant studies to rare variant studies for the same trait. Related to Figure 1 and Dataset EV1.

### `1C_DataCuration_Annotations.ipynb`
Curation of gene properties, and analysis of gene level overlap between common variant and rare variant associated genes for the same trait. Related to Figure 1 and Figure EV1. 

### `2_NetworkColocalization.ipynb`
Results of network colocalization analysis across 373 human traits. Related to Figure 2, and Dataset EV2. 

### `3_Reproducibility.ipynb`
Results of network colocalization analysis across 373 human traits, including results from multiple studies of the same trait, and analysis with multiple networks. Related to Figure 3, Figure EV2, Appendix Figure 3, and Dataset EV2. 

### `4A_FeaturePreparation.ipynb`
Curation of gene and trait features for each pair of common and rare variant studies. Related to Figure 4, Figure EV4, Dataset EV3 and Dataset EV4.

### `4B_FeatureAnalysis.ipynb`
Regression analysis of gene and trait features for prediction of network colocalization. Related to Figure 4, Figure 5, Figure EV3, Figure EV4, Appendix Figure 1, and Dataset EV4. 

### `5_Examples_NeuropsychiatricTraits.ipynb`
Identification and analysis of network colocalization for neuropsychiatric traits. Related to Figure 6 and Figure EV5.

### Supplemental Notebooks
**`S1_NetColocOptimization.ipynb`** - Optimization and benchmarking of NetColoc. Related to Appendix Figure 2 and Dataset EV5.   
**`S2_DilutionStudies.ipynb`** - Benchmarking of NetColoc in response to randomized inputs.  Related to Appendix Figure 3.
**`S3_ColocalizedNetwork_Uploads.ipynb`** - Formation of colocalized trait networks and upload to NDEx.    
**`S4_HCX_Creation.ipynb`** - Creation of hierarchical systems maps in HCX format.    

## CONTACT

For questions or issues, please open an issue on GitHub.