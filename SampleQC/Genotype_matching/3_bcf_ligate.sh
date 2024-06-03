#!/bin/bash

glimpse=/scratch/nsnyderm/programs/glimpse2
# full path to directory with GLIMPSE2 imputed bcfs per chunk
bcf_dir=$1
# full path to output directory for ligated bcf
ligated_bcf_outdir=$2
# final bcf output name in the format of ${output_prefix}.bcf 
ligated_bcf_outname=$3

ls -1v ${bcf_dir}/chr*.bcf > ${bcf_dir}/imputed_bcflist.txt

${glimpse}/GLIMPSE2_ligate_static --input ${bcf_dir}/imputed_bcflist.txt --output ${ligated_bcf_outdir}/${ligated_bcf_outname}
