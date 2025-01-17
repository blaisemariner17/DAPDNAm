#!/bin/bash

module load mamba/latest

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/nsnyderm/dap_rrbs/lid_pid`


# Array job : Allele-specific methylation

path_cgmap=/scratch/ccosta8/RRBS/Programs/xz_mod_cgmaptools/cgmaptools
path_ref=/scratch/bmarine2/Precision123-240506/genome/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa

bam_sort=/scratch/nsnyderm/dap_rrbs/sort_bam/${sampleID}.sort.bam

# extract sites to keep
module load vcftools-0.1.14-gcc-11.2.0
module load samtools-1.9-gcc-12.1.0
module load r-4.2.2-gcc-11.2.0

#index if the bam files do not already have an index .bai file
#samtools index $bam_sort

gunzip /scratch/bmarine2/Precision123-240506/QTL/sample_vcfs/${sampleID}.CGmap.PASS.DP5.vcf.gz

#get the cgmaptools commands that youve added to the bash rc file
source ~/.bashrc

# call allele-specific methylation/create asm files
cgmaptools asm -r $path_ref -b $bam_sort -l /scratch/bmarine2/Precision123-240506/QTL/sample_vcfs/${sampleID}.CGmap.PASS.DP5.vcf -m ass > /scratch/bmarine2/Precision123-240506/QTL/ASM/${sampleID}.CGmap.PASS.DP5.ASM
