// Note this does not curently support the --netfile input arguments
process getNetColocMatrices {

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

script:
def filter_arg = filter ? "--filter ${filter}" : ""
def alpha_arg  = alpha ? "--alpha ${alpha}" : ""

"""
python ${baseDir}/../carva/get_heat_matrix.py --outdir ${outdir} --uuid ${uuid} --name ${name} ${filter_arg} ${alpha_arg}

"""

}