#!/bin/bash

# Paths and directories
workingdir=/scratch/bmarine2/chrX-impute-241218
logdir=${workingdir}/logs/phase_bams
bamdir=/scratch/vsohrab/dap_map_canfam4/nvidia_fq2bam_mapping/bams
chunkdir=${workingdir}/glimpse_chunks
bcfoutdir=${workingdir}/glimpse_phase
mainscript=${workingdir}/4_impute_bams_glimps2_phase-BLM241219.sh

# Bamlist file
bamlist=${workingdir}/imputation_bamlist.txt

# Validate directories and files
if [ ! -d "$chunkdir" ]; then
    echo "Error: Chunk directory $chunkdir does not exist!"
    exit 1
fi

if [ ! -f "$bamlist" ]; then
    echo "Error: Bamlist file $bamlist does not exist!"
    exit 1
fi

# Determine the number of chunks dynamically
total_chunks=$(grep -w "chrX" ${chunkdir}/*chunks.chr*.txt | wc -l)
if [ "$total_chunks" -eq 0 ]; then
    echo "Error: No chunks found in $chunkdir for chrX!"
    exit 1
fi

# Submit SLURM array job
sbatch -p htc -q public --array=1-${total_chunks} \
    -t 0-04:00:00 --mem=120G --cpus-per-task 48 \
    -e ${logdir}/glimpse2_impute_phase_bams_gpu_%A_%a.e \
    -o ${logdir}/glimpse2_impute_phase_bams_gpu_%A_%a.o \
    ${mainscript} ${bcfoutdir} ${bamlist}

