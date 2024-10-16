## Quick guide to run files

----

`get_netcoloc_matrices.sh`
Calls *get_heat_matrix.py* to create saved versions of network heat matrix

----

`netcoloc_simulation.sh`
NOT AN SBATCH FILE
Runs a single trait or list of traits through netcoloc. 

----

`run_netcoloc_simulation.sh`
Runs a single netcoloc analysis via *do_netcoloc.py* where rare==common.

----

`overlap.sh`
Calculates overlap statistics between pairs of traits in parallel. Wrapper for *gene_overlap.py*

----

`run_create_sim_genesets.sh`
Wrapper for *create_sim_genesets.py*, creating gene sets with known features from a reference set

----

`run_netcoloc_between.sh`
Runs *do_netcoloc.py* for rare and common gene sets if rare != common.

----

`run_netcoloc_pairs.sh`
Runs *do_netcoloc.py* using a file of already defined pairs. Runs both RC and CR.

----

`run_netcoloc_within.sh`
Runs *do_netcoloc.py* for rare and common gene sets for a list of trait, with rare == common.

----

`run_overlap_pairs.sh`
Runs *gene_overlap.py* for trait pairs based on a file of already defined pairs. Runs both RC and CR. 

----

`run_shortest_paths.sh`
Runs *find_shortest_paths.py* with default options


