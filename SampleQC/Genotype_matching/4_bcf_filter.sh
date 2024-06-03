#!/bin/bash

set -e 

# load modules and list software
ml bcftools-1.14-gcc-11.2.0
plink2=/scratch/vsohrab/programs/plink2

raw_bcf=$1
filtered_bcf_outdir=$2

ref_fasta=/scratch/vsohrab/dap_map_canfam4/kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa
bcf_outname=`basename -s .bcf ${raw_bcf}`

## filter raw bcf file ##

# exclude variants/SNPs with max genotype probability < 0.7 in any sample within BCF
# split multi-allelic sites into biallelic entries
# exclude variants with incorrect or missing REF allele according to reference fasta file
# set variant ID to chr,pos,ref,alternate

bcftools view -Ou -e 'MAX(GP[*])<0.7' ${raw_bcf} |
        bcftools norm -m-any -Ou --threads 24 |
        bcftools norm -c wx -f ${ref_fasta} -Ou --threads 24 |
        bcftools annotate -Oz -I '%CHROM:%POS:%REF:%ALT' --threads 24 -o ${filtered_bcf_outdir}/${bcf_outname}.vcf.gz

echo "Finished filtering bcf file"

## create bcf index for filtered bcf file ##
bcftools index ${filtered_bcf_outdir}/${bcf_outname}.vcf.gz

echo "Finished creating index for filtered bcf file"

## create a plink set from filtered bcf file ##
${plink2} --vcf ${filtered_bcf_outdir}/${bcf_outname}.vcf.gz \
        --dog \
        --make-bed \
        --threads 24 \
        --out ${filtered_bcf_outdir}/${bcf_outname}_gp-0.70_biallelic

echo "Finished creating plink set from filtered bcf file"
