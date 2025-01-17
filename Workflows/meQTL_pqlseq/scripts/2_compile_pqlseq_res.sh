#!/bin/bash

awk 'FNR>1' pqlseq_results/pqlseq_res${SLURM_ARRAY_TASK_ID}_* > pqlseq_res_chr${SLURM_ARRAY_TASK_ID}_tmp.tsv
#cat header_pqlseq_res.tsv pqlseq_res_chr${SLURM_ARRAY_TASK_ID}_tmp.tsv > pqlseq_res_chr${SLURM_ARRAY_TASK_ID}.tsv
cat pqlseq_res_chr${SLURM_ARRAY_TASK_ID}_tmp.tsv > pqlseq_res_chr${SLURM_ARRAY_TASK_ID}.tsv
rm pqlseq_res_chr${SLURM_ARRAY_TASK_ID}_tmp.tsv
