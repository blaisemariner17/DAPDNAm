#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 48           # number of cores 
#SBATCH --mem=200G
#SBATCH -t 0-20:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)

module load bowtie2-2.4.2-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0
module load py-cutadapt-2.10-gcc-11.2.0

export PATH=$PATH:/scratch/nsnyderm/programs/FastQC:/scratch/nsnyderm/programs/bin:/scratch/nsnyderm/programs/TrimGalore-0.6.6:/scratch/nsnyderm/programs/bismark2/Bismark-0.24.0

sampleID=$1

ref_path=/scratch/bmarine2/Precision123-240506/genome/UU_Cfam_GSD_1.0-Y

out_path=/scratch/bmarine2/Precision123-240506/QTL/sample_bams
fastq_path=/scratch/nsnyderm/dap_rrbs/trimmed

cov_path=/scratch/bmarine2/Precision123-240506/QTL/cov_files

#bismark_genome_preparation --verbose ${ref_path}

bismark --genome_folder ${ref_path} --non_directional -o ${out_path} --score_min L,0,-0.6 -R 12 --parallel 24 ${fastq_path}/${sampleID}*trimmed.fq.gz

bismark_methylation_extractor -p -o ${cov_path} --no_overlap --bedGraph --comprehensive --parallel 10 --merge_non_CpG --genome_folder ${ref_path} ${out_path}/${sampleID}*.bam
