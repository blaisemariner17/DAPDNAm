#!/bin/bash

sampleID=LID_1069358_PID_10535

vcf_path=/scratch/bmarine2/Triad_mapping/nvidia_fq2bam_mapping/vcfs_merged/DAP_Triad_and_Precision_240602.noC.AF0.01.vcf.gz

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
	
sampleID=LID_1071411_PID_10559

vcf_path=/scratch/bmarine2/Triad_mapping/nvidia_fq2bam_mapping/vcfs_merged/DAP_Triad_and_Precision_240602.noC.AF0.01.vcf.gz

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
