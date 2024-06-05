#!/bin/bash
#SBATCH -G a100:2
#SBATCH -t 0-2
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12

sample=/scratch/bmarine2/Precision123-240506/QTL/lid_pids.txt


# Array job : Allele-specific methylation

path_cgmap=/scratch/ccosta8/RRBS/Programs/xz_mod_cgmaptools/cgmaptools
path_ref=/scratch/vsohrab/dap_map_canfam4/kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa
path_bam2=/scratch/nsnyderm/dap_rrbs/sort_bam/${sample}.sort.bam

# extract sites to keep
module load vcftools/0.1.12b
module load samtools/1.9

samtools index $path_bam2

gunzip /scratch/ccosta8/RRBS/Output/${sample}.R1_CGmap.PASS2.DP5.vcf.gz

# call allele-specific methylation/create asm files
$path_cgmap asm -r $path_ref -b $path_bam2 -l /scratch/bmarine2/Precision123-240506/QTL/sample_vcfs/${sample}.vcf -m ass > /scratch/bmarine2/Precision123-240506/QTL/ASM/${sample}.ASM
