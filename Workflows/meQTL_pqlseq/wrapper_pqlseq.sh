#!/bin/bash

for chr in {1..38};
do
#	if [[ "$chr" -eq 1 ]]
#	then
#		sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-527 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
#	fi
# 	if [[ "$chr" -eq 2 ]]
#       then
#           	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-544 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
#       fi
  	if [[ "$chr" -eq 3 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-355 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 4 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-336 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 5 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-508 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 6 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-331 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 7 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-241 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 8 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-271 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 9 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-373 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 10 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-396 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 11 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-253 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 12 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-207 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 13 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-223 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 14 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-149 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 15 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-203 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 16 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-299 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 17 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-243 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 18 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-347 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 19 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-137 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 20 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-276 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 21 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-187 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 22 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-141 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 23 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-154 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 24 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-369 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 25 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-309 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 26 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-314 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 27 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-151 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 28 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-349 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 29 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-129 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 30 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-206 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 31 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-206 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 32 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-61 --partition=htc -q public --time=0-00:30:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 33 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-107 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 34 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-217 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 35 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-130 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 36 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-88 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 37 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-128 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi
  	if [[ "$chr" -eq 38 ]]
        then
            	sbatch --job-name=pql_batch --mem=50G --cpus-per-task=20 --array=1-174 --partition=htc -q public --time=0-00:10:00 pqlseq_sbatch.sh $chr
        fi

done
