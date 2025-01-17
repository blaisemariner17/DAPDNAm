#!/bin/bash
sbatch --mem=50G --cpus-per-task=8 --array=1-1900 -t 0-1 2_sortbams.sh
