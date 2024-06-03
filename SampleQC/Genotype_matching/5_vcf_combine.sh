#!/bin/bash

set -e 

# load modules and list software
ml bcftools-1.14-gcc-11.2.0

bcftools merge $1 $2 -o DAP_Triad_and_Precision_240602.vcf.gz
