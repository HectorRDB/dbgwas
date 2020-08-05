# Running DBGWAS on docker

`docker run --rm leandroishilima/dbgwas:0.5.4 <args_to_DBGWAS>`

For mounting local folders to give the strains/databases as input and retrieve the output of DBGWAS, please use the `-v` docker parameter.

As an example, this is a commmand to run on the Pseudomonas Aeruginosa with Amikacin resistance phenotype:

`docker run --rm -v /data/leandro/docker_test/pseudomonas_aeruginosa_full_dataset:/pseudomonas_aeruginosa_full_dataset leandroishilima/dbgwas:0.5.4 -strains /pseudomonas_aeruginosa_full_dataset/strains -newick /pseudomonas_aeruginosa_full_dataset/strains.newick -nc-db /pseudomonas_aeruginosa_full_dataset/Resistance_DB_for_DBGWAS.fasta -pt-db /pseudomonas_aeruginosa_full_dataset/uniprot_sprot_bacteria_for_DBGWAS.fasta -output /pseudomonas_aeruginosa_full_dataset/output`


# Running DBGWAS on singularity

If you are having issues running DBGWAS on docker due to permissions, we recommend running it on singularity, e.g.:

`singularity run docker://leandroishilima/dbgwas:0.5.4 -strains pseudomonas_aeruginosa_full_dataset/strains -newick pseudomonas_aeruginosa_full_dataset/strains.newick -nc-db Resistance_DB_for_DBGWAS.fasta -pt-db uniprot_sprot_bacteria_for_DBGWAS.fasta` 

# Troubleshooting

`Nextflow` with `awsbatch` does not support images with entrypoints. In this case, you can use the no entrypoint image: `leandroishilima/dbgwas_no_entrypoint:0.5.4`.