#!/bin/bash
sbatch --mem=250G --cpus-per-task=8 --array=1 -t 0-4 4A_mergeASM.sh
