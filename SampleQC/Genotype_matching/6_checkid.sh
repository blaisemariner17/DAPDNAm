#!/bin/bash
#SBATCH -G a100:2
#SBATCH -t 0-12
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12

sampleIDs=/scratch/nsnyderm/dap_rrbs/lid_pid

vcf_path=/scratch/bmarine2/Triad_mapping/nvidia_fq2bam_mapping/DAP_Triad_and_Precision_240602.vcf.gz

#updated genome path for new dog genome
genome_path=/data/CEM/smacklab/genomes/bismark/UU_Cfam_GSD_1.0_ROSY/
fastq_trim=/scratch/nsnyderm/dap_rrbs/trimmed
sort_path=/scratch/nsnyderm/dap_rrbs/sort_bam

out_path=/scratch/bmarine2/Triad_mapping/mapped_all

module load bowtie2-2.4.2-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0
module load py-cutadapt-2.10-gcc-11.2.0

module load samtools-1.13-gcc-11.2.0

export PATH=$PATH:/scratch/nsnyderm/programs/FastQC:/scratch/nsnyderm/programs/bin:/scratch/nsnyderm/programs/TrimGalore-0.6.6:/scratch/nsnyderm/programs/bismark2/Bismark-0.24.0

for sampleID in $sampleIDs;
do
	
	if test -f ${sort_path}/${sampleID}.sort.bam
	then 
        	qtltools_path=/scratch/nsnyderm/programs/QTLtools_1.2_CentOS7.8_x86_64
        	${qtltools_path}/QTLtools mbv --vcf ${vcf_path} \
                 --bam ${sort_path}/${sampleID}.sort.bam \
                 --out check_out/${sampleID}.match
        exit
	fi

	echo "sorted bam does not exist, sorting bam first..."

	samtools sort -O BAM -@ 8 -o ${sort_path}/${sampleID}.sort.bam ${out_path}/${sampleID}*bam

	samtools index ${sort_path}/${sampleID}.sort.bam

	qtltools_path=/scratch/nsnyderm/programs/QTLtools_1.2_CentOS7.8_x86_64
	${qtltools_path}/QTLtools mbv --vcf ${vcf_path} \
        	--bam ${sort_path}/${sampleID}.sort.bam \
        	--out check_out/${sampleID}.match	
done




