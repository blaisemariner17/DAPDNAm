#!/bin/bash

sbatch --job-name=comp_pql --mem=50G --cpus-per-task=1 --array=10-38 --partition=htc -q public --time=0-00:30:00 2_compile_pqlseq_res.sh
