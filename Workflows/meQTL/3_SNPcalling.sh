#!/bin/bash

module load mamba/latest

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/nsnyderm/dap_rrbs/lid_pid`

#path_cgmap=/home/bmarine2/cgmaptools-0.1.3/bin/
bam_sort=/scratch/nsnyderm/dap_rrbs/sort_bam/${sampleID}.sort.bam

path_out=/scratch/bmarine2/Precision123-240506/QTL/sample_vcfs/${sampleID}.CGmap
path_ref=/scratch/bmarine2/Precision123-240506/genome/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa

#convert the bams to cgmap format
cgmaptools convert bam2cgmap -b $bam_sort -g $path_ref -o $path_out

#call SNPs from CGmap file
cgmaptools snv -i $path_out.ATCGmap.gz -m bayes -v $path_out.vcf --bayes-dynamicP -o $path_out.out --bayes-e=0.01 -a

#filter for DP>4 (read depth), remove "VAGUE" calls, and only include variants also present in Cayo genetic dataset (cayo.PASS.vcf.gz)
module load bcftools-1.10.2-gcc-12.1.0
#module load tabix/0.2.6

source activate snp_calling_env

bgzip -c $path_out.vcf > $path_out.vcf.gz
tabix -p vcf $path_out.vcf.gz
bcftools view -f 'PASS,.' $path_out.vcf.gz --output-type z > $path_out.PASS.vcf.gz
tabix -p vcf $path_out.PASS.vcf.gz
bcftools filter -i 'FORMAT/DP>4' $path_out.PASS.vcf.gz --output $path_out.PASS.DP5.vcf.gz --output-type z --regions-file /scratch/bmarine2/Triad_mapping/nvidia_fq2bam_mapping/vcfs_merged/DAP_Triad_and_Precision_240602.vcf.gz
tabix -p vcf $path_out.PASS.DP5.vcf.gz

rm $path_out.PASS.vcf.gz
rm $path_out.out
rm $path_out.ATCGmap.gz
