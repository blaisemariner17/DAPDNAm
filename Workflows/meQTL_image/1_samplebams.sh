#!/bin/bash
#SBATCH 


fastq_path=${1}
PID=${2}

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p ${fastq_path}/lids`

module load bowtie2-2.4.2-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0
module load py-cutadapt-2.10-gcc-11.2.0

export PATH=$PATH:/scratch/nsnyderm/programs/FastQC:/scratch/nsnyderm/programs/bin:/scratch/nsnyderm/programs/TrimGalore-0.6.6:/scratch/nsnyderm/programs/bismark2/Bismark-0.24.0

genome_path=/data/CEM/smacklab/genomes/bismark/UU_Cfam_GSD_1.0_ROSY/
out_path=/scratch/nsnyderm/dap_rrbs/mapped
cov_path=/scratch/nsnyderm/dap_rrbs/cov_files
trim_path=/scratch/nsnyderm/dap_rrbs/trimmed
finalcov_out=/scratch/nsnyderm/dap_rrbs/final_cov


[ -f ${cov_path}/${sampleID}_${PID}*.bismark.cov.gz ] && exit 0

## combine fastqs from mutliple lanes
cat ${fastq_path}/${sampleID}*R1*fastq.gz > /tmp/${sampleID}_${PID}.R1.fastq.gz

## trim adaptors
trim_galore --rrbs -j 8 \
        --non_directional \
        -o ${trim_path} --gzip \
        --basename ${sampleID}_${PID} \
        /tmp/${sampleID}_${PID}.R1.fastq.gz

rm /tmp/${sampleID}_${PID}.R1.fastq.gz

## map the reads
bismark --genome_folder ${genome_path} \
    -o ${out_path} \
    --score_min L,0,-0.6 -R 12 \
    --parallel 8 \
    --non_directional \
    ${trim_path}/${sampleID}_${PID}*.gz

## extract methylation data
bismark_methylation_extractor -s -o ${cov_path} \
        --gzip \
        --bedGraph --comprehensive --parallel 8 \
        --merge_non_CpG --genome_folder ${genome_path} \
        ${out_path}/${sampleID}_${PID}*.bam

##clean up unwanted files
rm ${cov_path}/*${sampleID}_${PID}*txt.gz
rm ${cov_path}/*${sampleID}_${PID}*bedGraph.gz

coverage2cytosine --merge_CpG --gzip --genome_folder $genome_path \
        -o ${sampleID}_${PID} \
        --dir ${finalcov_out}/ \
        ${cov_path}/${sampleID}_${PID}*.bismark.cov.gz
