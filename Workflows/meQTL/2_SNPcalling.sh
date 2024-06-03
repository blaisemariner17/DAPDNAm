#!/bin/bash
#SBATCH -G a100:2
#SBATCH -t 0-2
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12
path_cgmap=/scratch/ccosta8/RRBS/Programs/xz_mod_cgmaptools/cgmaptools

sample=/scratch/bmarine2/Precision123-240506/QTL/lid_pids.txt

path_bam=/scratch/bmarine2/Precision123--240506/QTL/sample_bams/${sampleID}*.bam
path_bam_sort=/scratch/bmarine2/Precision123--240506/QTL/sample_bams/sort/${sampleID}*.bam
path_out=/scratch/bmarine2/Precision123--240506/QTL/sample_vcfs/${sampleID}*.CGmap
path_ref=/scratch/vsohrab/dap_map_canfam4/kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa

#sort the bam files
samtools sort -o $path_bam_sort $path_bam
samtools view $path_bam | wc -l

#convert the bams to cgmap format
$path_cgmap convert bam2cgmap -b $path_bam_sort -g $path_ref -o $path_out

#filter for DP>4 (read depth), remove "VAGUE" calls, and only include variants also present in Cayo genetic dataset (cayo.PASS.vcf.gz)
module bcftools/1.10.2
module load tabix/0.2.6

bgzip -c $path_out.vcf > $path_out.vcf.gz
tabix -p vcf $path_out.vcf.gz
bcftools view -f 'PASS,.' $path_out.vcf.gz --output-type z > $path_out.PASS.vcf.gz
tabix -p vcf $path_out.PASS.vcf.gz
bcftools filter -i 'FORMAT/DP>4' $path_out.PASS.vcf.gz --output $path_out.PASS.DP5.vcf.gz --output-type z --regions-file /scratch/nsnyderm/meQTL/cayo.PASS.vcf.gz
tabix -p vcf $path_out.PASS.DP5.vcf.gz

rm $path_out.PASS.vcf.gz
rm $path_out.out
rm $path_out.ATCGmap.gz
