#!/bin/bash
module load bcftools-1.14-gcc-11.2.0
module load htslib-1.16-gcc-11.2.0
module load samtools-1.16-gcc-11.2.0

## Argument 1 is the full path for the output of imputed BCFs
## Argument 2 is a list of BAM file paths and the animal ID (space delimited)

set -e

## Prep for GLIMPSE2
path_to_glimpse=/scratch/nsnyderm/programs/glimpse2
path_to_dog10k_imputation_ref=/scratch/vsohrab/dap_glimpse_impute_dog10k
path_to_glimpse2_chunks=/scratch/bmarine2/chrX-impute-241218/glimpse_chunks
path_to_glimpse2_binary_ref=${path_to_glimpse2_chunks}/glimpse_binary_ref_panel
imputed_bcf_outdir=${1}
bamlist=${2}

# Chromosome
chr="chrX"

# Get the line for the current task ID
LINE=$(grep -w "chrX" ${path_to_glimpse2_chunks}/*.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Safety check for empty line
if [ -z "$LINE" ]; then
    echo "Error: No chunk found for task ID ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# Extract each chunk entry information
CHR=$(echo $LINE | cut -d" " -f2)   # Chromosome and range
IRG=$(echo $LINE | cut -d" " -f3)  # Start position of chunk
ORG=$(echo $LINE | cut -d" " -f4)  # End position of chunk
REGS=$IRG                          # Start of the range
REGE=$ORG                          # End of the range

# Debugging output
echo "Chunk info: $LINE"
echo "CHR=$CHR, IRG=$IRG, ORG=$ORG, REGS=$REGS, REGE=$REGE"

# Check if final output already exists
#test -f ${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}.bcf.csi && exit

# Filter male and female BAM files
# Filter male and female BAM files
# Filter male and female BAM files
male_bamlist=/scratch/bmarine2/chrX-impute-241218/male_bamlist.txt
female_bamlist=/scratch/bmarine2/chrX-impute-241218/female_bamlist.txt

# Debugging output for BAM lists
#echo "Male BAM list:"
#cat ${male_bamlist}
#echo "Female BAM list:"
#cat ${female_bamlist}

# Start imputation for male samples
male_bcf=${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}_male.bcf
if [ -s ${male_bamlist} ]; then
    echo "Starting imputation for male samples at $(date '+%Y-%m-%d %H:%M:%S')"
    ${path_to_glimpse}/GLIMPSE2_phase_static --bam-list ${male_bamlist} \
        --reference /scratch/bmarine2/chrX-impute-241218/glimpse_binary_ref_panel/Dog10K_chrX_1_6603974.bin \
        --threads 48 \
        --output ${male_bcf}
    echo "Imputation for male samples has completed at $(date '+%Y-%m-%d %H:%M:%S')"
else
    echo "No male samples found for this chunk. Skipping."
    male_bcf=""
fi

# Start imputation for female samples
female_bcf=${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}_female.bcf
if [ -s ${female_bamlist} ]; then
    echo "Starting imputation for female samples at $(date '+%Y-%m-%d %H:%M:%S')"
    ${path_to_glimpse}/GLIMPSE2_phase_static --bam-list ${female_bamlist} \
        --reference /scratch/bmarine2/chrX-impute-241218/glimpse_binary_ref_panel/Dog10K_chrX_1_6603974.bin \
        --threads 48 \
        --output ${female_bcf}
    echo "Imputation for female samples has completed at $(date '+%Y-%m-%d %H:%M:%S')"
else
    echo "No female samples found for this chunk. Skipping."
    female_bcf=""
fi

# Combine male and female BCFs
if [ -n "$male_bcf" ] && [ -n "$female_bcf" ]; then
    echo "Combining male and female BCFs..."
    bcftools merge -O b -o ${imputed_bcf_outdir}/male_female_combined_chrX.bcf ${male_bcf} ${female_bcf}
    echo "Combined BCF created: ${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}.bcf"
elif [ -n "$male_bcf" ]; then
    echo "No female samples found. Using male BCF as the output."
    mv ${male_bcf} ${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}.bcf
elif [ -n "$female_bcf" ]; then
    echo "No male samples found. Using female BCF as the output."
    mv ${female_bcf} ${imputed_bcf_outdir}/${CHR}_${REGS}_${REGE}.bcf
else
    echo "Error: No male or female BCFs generated for chunk $CHR:${REGS}-${REGE}!"
    exit 1
fi

bcftools index -f ${imputed_bcf_outdir}/male_female_combined_chrX.bcf

# Clean up temporary files
rm -f ${male_bamlist} ${female_bamlist}

