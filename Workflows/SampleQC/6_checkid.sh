#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 48           # number of cores 
#SBATCH --mem=200G
#SBATCH -t 0-20:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)

sampleIDs_=/scratch/bmarine2/Precision123-240506/QTL/lid_pids.txt
sampleID=(`awk '{print $1}' $sampleIDs_`)

vcf_path=/scratch/bmarine2/Triad_mapping/nvidia_fq2bam_mapping/DAP_Triad_and_Precision_240602.vcf.gz

#updated genome path for new dog genome
fastq_trim=/scratch/nsnyderm/dap_rrbs/trimmed
sort_path=/scratch/nsnyderm/dap_rrbs/sort_bam

out_path=/scratch/bmarine2/Triad_mapping/checkid

module load bowtie2-2.4.2-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0
module load py-cutadapt-2.10-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0

export PATH=$PATH:/scratch/nsnyderm/programs/FastQC:/scratch/nsnyderm/programs/bin:/scratch/nsnyderm/programs/TrimGalore-0.6.6:/scratch/nsnyderm/programs/bismark2/Bismark-0.24.0

qtltools_path=/scratch/nsnyderm/programs/QTLtools_1.2_CentOS7.8_x86_64
	${qtltools_path}/QTLtools mbv --vcf ${vcf_path} \
	--bam ${sort_path}/${sampleID}.sort.bam \
	--out ${out_path}/${sampleID}.match
