#!/bin/bash

# Directories
workingdir=/scratch/bmarine2/chrX-impute-241218
bamdir=/scratch/vsohrab/dap_map_canfam4/nvidia_fq2bam_mapping/bams
bamlist=${workingdir}/imputation_bamlist.txt
tsv_file=/scratch/bmarine2/chrX-impute-241218/all_dogs_mappedto_Sex.tsv

# Create an empty bamlist
> $bamlist

# Loop through each ID in the first column of the TSV file
while IFS=$'\t' read -r animal_id sex; do
    # Find the BAM file that starts with the animal ID (matching the structure in bamdir)
    bam_file=$(find $bamdir -name "${animal_id}.canfam4.bam" -print -quit)

    # If a BAM file was found, add the full path and ID to the bamlist
    if [ -n "$bam_file" ]; then
        echo "$bam_file $animal_id" >> $bamlist
    fi
done < $tsv_file

echo "Bamlist has been created: $bamlist"

