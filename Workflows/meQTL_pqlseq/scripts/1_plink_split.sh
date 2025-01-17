#!/bin/bash

# Load necessary modules
module load mamba/latest
module load plink/1.90
module load bedtools2-2.30.0-gcc-11.2.0
module load samtools-1.9-gcc-12.1.0
module load bcftools-1.14-gcc-11.2.0
module load tabix-2013-12-16-gcc-12.1.0
module load htslib-1.16-gcc-11.2.0

chr=chr${SLURM_ARRAY_TASK_ID}

# Set paths to files and directories
methylation_matrix="/scratch/bmarine2/meQTL-241101-GEMMA/perc_meth_${chr}.txt"  # Percent methylation matrix
methylation_matrix_forgemma="/scratch/bmarine2/meQTL-241101-GEMMA/perc_meth_forgemma_${chr}.txt"
output_dir="/scratch/bmarine2/meQTL-241101-GEMMA/gemma_out"                   # Directory for output files
gemma_path="/scratch/bmarine2/meQTL-241101-GEMMA/gemma-0.98.5-linux-static-AMD64"
bcf_imputed="/scratch/vsohrab/dap_glimpse_impute_dog10k/glimpse_ligate/DogAgingProject_2023_N-7627_canfam4.bcf"
plink_base="${output_dir}/filtered_genotypes_${chr}"                           # Base name for PLINK files
vcf_file="${output_dir}/combined_${chr}.vcf"                                   # Combined VCF for chr14

echo $chr
echo "Output directory is ${output_dir}"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

echo 'Generate a XX kb window around each methylation site'
awk -F"_" '{start=$2-50000; end=$3+50000; if (start < 0) start=1; print $1"\t"start"\t"end}' "$methylation_matrix" > "${output_dir}/methylation_regions_${chr}.bed"

echo 'Get the vcf file from the imputed vcf'
#sample_list="/scratch/bmarine2/meQTL-241101-GEMMA/methylation_samples.txt"
#bcftools view ${bcf_imputed} --regions $chr --samples-file $sample_list --threads 4 > ${vcf_file}

echo 'Filter VCF for SNPs within 50kb of CpG regions'
bedtools intersect -b "$vcf_file" -a "${output_dir}/methylation_regions_${chr}.bed" > "${output_dir}/filtered_regions_${chr}.vcf"
#bedtools intersect -a "$vcf_file" -b "${output_dir}/methylation_regions_${chr}.bed" > "${output_dir}/filtered_${chr}.vcf"


## add the header back after the bcftools intersect call
grep "^#" "$vcf_file" > "${output_dir}/header_with_samples.txt"
cat "${output_dir}/header_with_samples.txt" "${output_dir}/filtered_${chr}.vcf" > "${output_dir}/filtered_${chr}_header.vcf"

echo 'Convert the filtered VCF to PLINK format'
plink --vcf "${output_dir}/filtered_${chr}_header.vcf" --make-bed --out "$plink_base" --allow-extra-chr --chr-set 38 --geno 0 --snps-only --maf 0.1 --r2 --ld-window-kb 70 --ld-window-r2 0.2
