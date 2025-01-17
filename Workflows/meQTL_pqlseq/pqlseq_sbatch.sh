#!/bin/bash

module load r-4.2.2-gcc-11.2.0
Rscript meqtl_pqlseq-241101.R $1 $SLURM_ARRAY_TASK_ID

##can also run this in 2_xxxx.sh file
#awk 'FNR>1' pqlseq_results/pqlseq_res${SLURM_ARRAY_TASK_ID}_* > pqlseq_res_chr${SLURM_ARRAY_TASK_ID}_tmp.tsv
#cat pqlseq_res_chr${SLURM_ARRAY_TASK_ID}_tmp.tsv > pqlseq_res_chr${SLURM_ARRAY_TASK_ID}.tsv
#rm pqlseq_res_chr${SLURM_ARRAY_TASK_ID}_tmp.tsv
