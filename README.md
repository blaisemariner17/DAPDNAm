# DAPDNAm
A GitHub R package with functions and workflows for the processing of dog <ins>D</ins>NA <ins>m</ins>ethylation data for the <ins>D</ins>og <ins>A</ins>ging <ins>P</ins>roject (DAP). More information on DAP can be found at <ins>dogagingproject.org</ins>.

### Installation

```
devtools::install_github("blaisemariner17/DAPDNAm")
```

## Directories of interest

Workflows contain the workflows that one can use to analyze their own data. Differetial_analysis includes all of the sample metadata creation, region metadata creation, region exploration, differential methylation with a binomial mixed modeling approach, and over-representation analysis. SampleQC contains the scripts that take your fq files to sample matching your bismark files (checkid.sh, being the final step). Then, you can use the provided R script (8_dap_check_lid.R) to see if there are any genotyping mismatches. The meQTL directory provides the scripts needed to do methylation QTL on your samples in numbered order.

The functions that are provided in the R function (and available for use after installing DAPDNAm on your machine), are used in the differential analysis workflows.
