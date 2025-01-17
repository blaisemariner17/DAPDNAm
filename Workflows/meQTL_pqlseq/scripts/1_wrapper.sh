#!/bin/bash

sbatch --job-name=plink_split --mem=50G --cpus-per-task=4 --array=1-38 --partition=htc -q public --time=0-02:00:00 1_plink_split.sh
