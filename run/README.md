# Guide to run files

## Gene Set Overlaps

`rvc_overlap.sh` calculates the direct overlap between genes identified by rare variants and common variants for each trait.

* Inputs:
    * `<rare_file>` - file listing trait names/file prefixes of rare variant gene sets (one per line)
    * `<common_file>` - file listing trait names/file prefixes of common variant gene sets (one per line)
    * (Optional) P-value thresholds, gene background size, minimum number of genes per set
* Executes: `gene_overlap.py`
* Outputs: Results of hypergeometric tests for each trait pair of rare and common gene sets, written to OUTDIR. 

## Pre-calculation of heat matrices for each network

Calculation of the heat matrix for a network is independent of the gene sets being analyzed. Precomputing these matrices improves efficiency of NetColoc across many gene sets. The script `get_netcoloc_matrices.sh` performs this step for multiple networks.

* Inputs: Network name (for file naming) and network UUID from [NDEx](https://ndexbio.org)
* Executes: `get_heat_matrix.py`
* Outputs:
    * `<name>_w_prime.npy` - Normalized adjacency matrix
    * `<name>_individual_heats.npy` - Heat diffusion matrix for network propagation
    * `<name>_nodes.txt` - List of all network nodes
    * `<name>_degrees.txt` - Node degrees

## Network Colocalization

`rvc_qnetcoloc.sh` performs the main network colocalization analysis.

* Inputs:
    * `<config>` - Configuration file for NetColoc (see `run_configs/example.config`) specifying data directories, trait lists, and NetColoc options
    * Network information: `<uuid>` network identifier, `<name>` network name to use for file naming (must match pre-calculated matrix files), and `<netdir>` directory containing precomputed heat matrix
    * Gene list suffixes: `<raresuff>`, `<commonsuff>` for file naming
    * Run parameters: `<binsize>` for degree-matched randomization, `<min_genes>` minimum genes per set
* Executes: `do_carva_netcoloc.py`
* Outputs:
    * Z-scores: `<trait>_z_<RV|CV>_q_<transform>_<norm>.tsv` one file per gene set
    * NetColoc results: `[q]netcoloc_<traitR>_<traitC>_<suffix>.txt` one file per trait pair

## Network Shortest Paths

`run_shortest_paths.sh` calculates shortest paths between all network nodes.

* Inputs: Network UUID, output prefix for file naming
* Executes: `find_shortest_paths.py`
* Outputs: `<prefix>path_lengths.csv` containing all-pairs shortest path lengths

## Network Gene Set Properties

`annotations.sh` calculates network topology properties of gene sets. Requires shortest paths to be previously calculated using `run_shortest_path.sh`.

* Inputs: 
    * `<geneset_file>` - file listing file names of rare variant gene sets 
    * `<geneset_file2>` - file listing file names of common variant gene sets
    * `<uuid_list>` - list of network identifiers on NDEx
    * `<name_list>` - names for each network in uuid_list
* Executes: 
    * `get_network_stats.py` - calculate all clustering coefficients for each node
    * `network_annotation.py` - calculate all network annotations
* Outputs:
    * `network_stats_...` - Network statistics for all genesets (except modularities)
    * `network_modularities...` - Modularities for all genesets

## Subnetwork Creation

`run_subnetwork_creation.sh` generates subnetworks from network colocalization results.

* Inputs:
    * `<traitlist_file>` - TSV file with `TraitR` and `TraitC` columns listing trait pairs
    * `<z_dir>` - Directory containing z-scores from network propagation
    * `<uuid>` - Network UUID for subnetwork extraction
    * (Optional) Z-score thresholds for defining colocalized genes
* Executes: `create_subnetworks.py`
* Outputs:
    * `<traitR>_<traitC>_subnetwork_all.tsv` - Edgelist with input + colocalized genes
    * `<traitR>_<traitC>_subnetwork_all_node_attributes.tsv` - Node attributes
    * `<traitR>_<traitC>_subnetwork_inputs.tsv` - Edgelist with input genes only
    * `<traitR>_<traitC>_subnetwork_inputs_node_attributes.tsv` - Node attributes 

## Download Summary Statistics

The files `download_sumstats.sh` and `run_download_sumstats.sh` can be used for bulk downloads of summary statistics from the GWAS Catalog.  

* Inputs: File containing download urls (e.g. http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005842)
* Outputs: Summary statistics files

## Additional Analysis Scripts

**Network Shuffling**: `network_shuffle.sh` and `rvc_shuffled_networks.sh` generate degree-preserving randomized networks for null model comparisons. See `network_shuffle.py` documentation for details.

**Simulated Gene Sets**: `run_create_sim_genesets.sh` and `run_create_sim_genesets_quant.sh` create simulated gene sets with controlled overlap and relevance for benchmarking. Execute `create_sim_genesets.py` with specified parameters.

**Alternative Networks**: `rvc_qnetcoloc_othernet.sh` runs NetColoc analysis using networks other than PCNet 2.0.

**Parameter Optimization**: `sweep_th_rvc_qnetcoloc.sh` and `sweep_bin_rvc_qnetcoloc.sh` perform parameter optimization of Q-NetColoc by executing `do_carva_netcoloc.py` over a range of threshold and binning parameters.



