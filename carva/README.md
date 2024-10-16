## Quick guide to modules/executables

### Modules

`data_utils.py` - functions for working with meta data
`geneset_utils.py` - functions for loading gene sets
`network_utils.py` - functions for loading and working with networks

### Executables
#### `create_sim_genesets.py`

**Inputs:**
* file of original gene sets
* file of network nodes
* desired overlap
* desired relevance
* desired total number of genes
* number of repeats to create
* background

**Computation:**
* Filters geneset based on node presence in network.
* Filters to genesets with >= desired number of total genes
* Partition gene set into overlap, unique1 and unique2
* If relevance < 1, add noise to the sets (DOES NOT add noise to overlap genes). If background is degree, add noise nodes with similar degree. 

**Outputs:**
* For each repeat of each geneset, a simulated set of 'CV' and 'RV' genes, written to file  

----

#### `do_netcoloc.py`

**Inputs:**
* I/O directories
* rare and common trait names
* network information
* method for overlap control ('remove' or 'bin')

**Computation:**
* Load gene sets and checks for at least 3 genes in each set (before and after filtering to network)
* Subset to top 500 if > 500 genes in set
* Check if z-scores already exist, otherwise calculate using `netcoloc`
* write all netcoloc metrics and permutation test metrics to file. 

**Outputs:**

------
#### `find_shortest_paths.py`

Loads a network from NDEx or file, and calculates all shortest paths between nodes in that network.

----

#### `gene_overlap.py`

**Inputs:**
* rare and common trait names
* data paths
* rare and common thresholds
* minimum genes
* background N
* test name

**Computation:**
Identifies shared genes between the two sets and perform a hypergeometric test

**Outputs:**
Statistics written to file

----

#### `get_heat_matrix.py`

Calls `netcoloc` functions to create and save `w_prime` and the `individual_heats` for easy use across many tests

