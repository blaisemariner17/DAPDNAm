#!/bin/bash
#SBATCH -G a100:2
#SBATCH -t 0-2
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12

module load anaconda/py3
source activate bismark
module load trim_galore/0.4.0

sampleID=/scratch/bmarine2/Precision123-240506/QTL/lid_pids.txt
genome_path=/scratch/vsohrab/dap_map_canfam4/kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa

out_path=/scratch/bmarine2/Precision123--240506/QTL/sample_bams
fastq_path=/scratch/nsnyderm/dap_rrbs/trimmed

bismark --genome_folder ${genome_path} --non_directional -o ${out_path} --score_min L,0,-0.6 -R 12 --parallel 4 ${fastq_path}/${sampleID}*trimmed.fq.gz

bismark_methylation_extractor -p -o ${cov_path} --no_overlap --bedGraph --comprehensive --parallel 10 --merge_non_CpG --genome_folder ${genome_path} ${out_path}/${sampleID}*.bam

conda deactivate
