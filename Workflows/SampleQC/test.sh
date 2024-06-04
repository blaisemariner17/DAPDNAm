#!/bin/bash

sampleIDs_=/scratch/bmarine2/Precision123-240506/QTL/lid_pids.txt
sampleIDs=(`awk '{print $1}' $sampleIDs_`)

for i in "${!sampleIDs[@]}";
do
	echo "${sampleIDs[i]}"
done
