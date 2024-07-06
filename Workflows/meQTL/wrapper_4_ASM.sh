#!/bin/bash
sbatch --mem=50G --cpus-per-task=8 --array=1-5 -t 0-1 4_ASM_RRBS.sh
