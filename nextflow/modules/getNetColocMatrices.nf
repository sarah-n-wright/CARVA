// Note this does not curently support the --netfile input arguments
process getNetColocMatrices {

errorStrategy { task.attempt <= task.maxRetries ? 'retry' : 'ignore' }
maxRetries 3
conda 'netcoloc.yml'

input: 
    val uuid
    val outdir
    val name
    val filter
    val alpha

output:
    path "${name}_w_prime.npy", emit: w_prime
    path "${name}_individual_heats.npy", emit: individual_heats
    path "${name}_nodes.txt", emit: nodes
    path "${name}_degrees.txt", emit: degrees
    path "netcoloc_matrices", emit: netdir

script:
def filter_arg = filter ? "--filter ${filter}" : ""
def alpha_arg  = alpha ? "--alpha ${alpha}" : ""

"""
python ${baseDir}/../carva/get_heat_matrix.py --outdir ${outdir} --uuid ${uuid} --name ${name} ${filter_arg} ${alpha_arg}
mkdir netcoloc_matrices
cp ${name}_w_prime.npy ${name}_individual_heats.npy ${name}_nodes.txt ${name}_degrees.txt netcoloc_matrices/
"""

}