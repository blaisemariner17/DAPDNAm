#!/bin/bash

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/nsnyderm/dap_rrbs/lid_pid`

module load samtools/1.9

bams=/scratch/nsnyderm/dap_rrbs/mapped/${sampleID}_trimmed_bismark_bt2.bam
bam_sort=/scratch/nsnyderm/dap_rrbs/sort_bam/${sampleID}.sort.bam

#sort the bam files
samtools sort -o $bam_sort $bams
#samtools view $path_bam | wc -l
