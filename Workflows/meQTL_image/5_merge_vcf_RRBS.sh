#!/bin/sh

module load bcftools-1.10.2-gcc-12.1.0

bcftools merge /scratch/bmarine2/Precision123-240506/QTL/sample_vcfs/*.CGmap.PASS.DP5.vcf.gz --output /scratch/bmarine2/Precision123-240506/QTL/sample_vcfs/all_RRBS_samples.CGmap.PASS.DP5.vcf.gz --output-type z --force-samples
