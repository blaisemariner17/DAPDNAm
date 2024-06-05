#!/bin/bash

sampleIDs=(`awk '{print $1}' /scratch/bmarine2/Precision123-240506/QTL/lid_pids.txt`)
#echo ${sampleIDs//\"/}
sbatch 1_sampleVCFs.sh ${sampleIDs//\"/}
