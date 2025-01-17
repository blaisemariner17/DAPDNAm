#!/bin/bash

#SBATCH -N 1            # number of nodes
#SBATCH -c 5           # number of cores 
#SBATCH --mem=200G
#SBATCH -t 0-20:00:00   # time in d-hh:mm:ss
#SBATCH -p general      # partition 
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --mail-user="bmarine2@asu.edu"
#SBATCH --export=NONE   # Purge the job-submitting shell environment


module load r-4.2.2-gcc-11.2.0
Rscript meQTL.R
