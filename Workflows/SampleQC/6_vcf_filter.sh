#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 30          # number of cores 
#SBATCH --mem=200G
#SBATCH -t 0-3:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)

# load modules and list software
ml bcftools-1.14-gcc-11.2.0

bcftools view -m2 -M2 -v snps -e 'ALT="C" || REF="C" || INFO/AF<0.01' --threads 24 -Oz -o /scratch/bmarine2/Triad_mapping/nvidia_fq2bam_mapping/vcfs_merged/DAP_Triad_and_Precision_240602.noC.AF0.01.vcf.gz /scratch/bmarine2/Triad_mapping/nvidia_fq2bam_mapping/vcfs_merged/DAP_Triad_and_Precision_240602.vcf.gz
