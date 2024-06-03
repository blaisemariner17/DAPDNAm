#!/bin/bash  
                                                                 
module load bcftools-1.14-gcc-11.2.0
module load htslib-1.16-gcc-11.2.0

## argument 1 is the full path for the output of imputed bcfs
## argument 2 is a list of bam files paths and the animal ID (space delimited)

set -e

## prep for GLIMPSE2
path_to_glimpse=/scratch/nsnyderm/programs/glimpse2
path_to_dog10k_imputation_ref=/scratch/vsohrab/dap_glimpse_impute_dog10k
path_to_glimpse2_chunks=${path_to_dog10k_imputation_ref}/glimpse_chunks
path_to_glimpse2_binary_ref=${path_to_dog10k_imputation_ref}/glimpse_binary_ref_panel
imputed_bcf_outdir=${1}
bamlist=${2}
LINE=`cat ${path_to_glimpse2_chunks}/*chunks.chr*.txt | sed -n ${SLURM_ARRAY_TASK_ID}p` 

# extract each chunk entry information
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
CHR=$(echo ${LINE} | cut -d" " -f2)
REGS=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f1)
REGE=$(echo ${IRG} | cut -d":" -f 2 | cut -d"-" -f2)
REF="Dog10K" # reference binary panel prefix

# if bcf for this chunk already exists, exit script
test -f ${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}.bcf.csi && exit


echo "Starting imputation for $LINE at $(date '+%Y-%m-%d %H:%M:%S')"
 
# impute
${path_to_glimpse}/GLIMPSE2_phase_static --bam-list ${bamlist} \
        --reference ${path_to_glimpse2_binary_ref}/${REF}_${CHR}_${REGS}_${REGE}.bin \
        --threads 48 \
        --output ${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}.bcf

echo "Imputation for $LINE has completed at $(date '+%Y-%m-%d %H:%M:%S')"
