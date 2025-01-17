#!/bin/bash
module load samtools-1.16-gcc-11.2.0
module load bcftools-1.14-gcc-11.2.0
module load htslib-1.16-gcc-11.2.0

# Define paths
path_to_glimpse=/scratch/nsnyderm/programs/glimpse2
path_to_glimpse2_binary_ref=/scratch/bmarine2/chrX-impute-241218/glimpse_binary_ref_panel
path_bcf=/scratch/vsohrab/dap_glimpse_impute_dog10k
vcf_ref_panel=${path_to_glimpse2_binary_ref}/AutoAndXPAR.Dog10K.chrX.phased.sites.vcf.gz
path_to_glimpse_chunks=/scratch/bmarine2/chrX-impute-241218/glimpse_chunks

chrom=X

# Step 1: Extract chrX reference VCF
bcftools view -G -Oz -r chr${chrom} -o ${vcf_ref_panel} \
              ${path_bcf}/AutoAndXPAR.Dog10K.Phased.AddAlleleNumber.bcf
bcftools index -f ${vcf_ref_panel}

# Step 2: Create GLIMPSE chunks
${path_to_glimpse}/GLIMPSE2_chunk_static \
    --input ${vcf_ref_panel} \
    --region chr${chrom} --output ${path_to_glimpse_chunks}/Dog10K.chunks.chr${chrom}.txt --window-size 5000000 --overlap-size 100000

# Check if chunk file exists and is non-empty
if [ ! -s ${path_to_glimpse_chunks}/Dog10K.chunks.chr${chrom}.txt ]; then
    echo "Error: Chunk file is empty or missing." >&2
    exit 1
fi

echo "Finished chunking for chr${chrom} at $(date '+%Y-%m-%d %H:%M:%S')"

# Step 3: Prepare GLIMPSE binary reference panel
REF=${path_bcf}/AutoAndXPAR.Dog10K.Phased.AddAlleleNumber.bcf
while IFS="" read -r LINE || [ -n "$LINE" ];
do
    printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
    IRG=$(echo $LINE | cut -d" " -f3)
    ORG=$(echo $LINE | cut -d" " -f4)

    ${path_to_glimpse}/GLIMPSE2_split_reference_static --reference ${REF} \
        --input-region ${IRG} --output-region ${ORG} \
        --output ${path_to_glimpse2_binary_ref}/Dog10K
done < ${path_to_glimpse_chunks}/Dog10K.chunks.chr${chrom}.txt

awk '{print $2, $3, $6, $7}' OFS='\t' glimpse_chunks/Dog10K.chunks.chrX.txt > glimpse_chunks/Dog10K.chunks.chrX.fixed.txt

echo "Finished creating GLIMPSE2 binary reference panel for chr${chrom} at $(date '+%Y-%m-%d %H:%M:%S')"

