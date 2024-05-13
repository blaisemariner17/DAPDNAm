#### Blaise Mariner now has a slight grasp on Sol and the dap mDNA data so now it is time to make a full workflow pipeline
#### This is part 1 of this series, which prepared the regions, their coverage and methylation, 
#### so that we can visualize them in pt 2
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2023-11-10
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project
rm(list=ls())

#### GLOBALS ####
# this scrip errored with 400gb 12 cores parallelized
# takes 18 minutes with 6 cores
cores_ <- 8

#### END GLOBALS ####

library_list <- c(
  "ggplot2",
  "svglite",
  "ggExtra",
  "ggtext",
  "GenomicRanges",
  "GenomicFeatures",
  "tidyverse",
  "bsseq",
  "comethyl",
  "PQLseq",
  "foreach",
  "parallel",
  "Biostrings"
)

lapply(library_list, require, character.only = TRUE)

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

## metadata
#just the precision dogs used here
p1 <- read_rds("metadata_samples/P1-DAP-metaData-240510.rds")
p2 <- read_rds("metadata_samples/P2-DAP-metaData-240510.rds")
p3 <- read_rds("metadata_samples/P3-DAP-metaData-240510.rds")

metaData <- rbind(p1,p2,p3);nrow(metaData)

#### do the following, or skip to the "OR:" section to save time ####
#load in the rds file
dap <- readRDS("/scratch/nsnyderm/dap_rrbs/bismarkBSseq.rds")
dap
sampleNames(dap) <- paste(sub(".*.(LID_\\d+)_.*", "\\1", sampleNames(dap)), 
                          sub(".*_(PID_\\d+).*", "\\1", sampleNames(dap)), 
                          sep = "_")

dap <- dap[,sampleNames(dap) %in% rownames(metaData)]

metaData <- metaData[metaData$lid_pid %in% sampleNames(dap),]

## split up by chromosome into a list
chrs=paste0("chr",c(1:38,"X"))

dap_list=parallel::mclapply(chrs,
                            function(x){
                              chrSelectBSseq(dap, seqnames = x, order = TRUE)
                            },
                            mc.cores = 8)

names(dap_list)=paste0("chr",c(1:38,"X"))

cpg_filtering <- function(x, metaData){
  x_oi <- x[,sampleNames(x) %in% metaData$lid_pid[metaData$Cohort == "precision_1"]]
  getcov <- getCoverage(x_oi, type = "Cov")
  keep = rowSums(getcov>=5)>=(.33*ncol(x_oi))
  return(x[keep,])
}

dap_list_filtered=parallel::mclapply(dap_list,   
                                     cpg_filtering,
                                     metaData,
                                     mc.cores = 12)

names(dap_list_filtered)

res <- 0
for(i in 1:length(dap_list_filtered)){
  res <- res + length(dap_list_filtered[[i]])
}
res
rm(dap)

if (dir.exists("getcoverage_data_and_plots") == FALSE){
  dir.create("getcoverage_data_and_plots")
}

# devtools::install_github("blaisemariner17/DAPDNAm", force = TRUE)
library("DAPDNAm")
# DAPDNAm::filtering_region_and_coverage(dap_list_filtered[['chr1']])

for (maxgap in c(250)){
  print(maxgap)
  time_start <- Sys.time()
  print(time_start)
  #remotes::install_version("matrixStats", version="1.1.0") #use this if you get a weird error in the getCoverage that is unclear
  dap_reduced=parallel::mclapply(dap_list_filtered,
                                 DAPDNAm::filtering_region_and_coverage,
                                 metaData = metaData,
                                 maxgap = maxgap,
                                 in_at_least_X_samples = .33*ncol(dap_list_filtered[[1]]),
                                 at_least_X_coverage_of_a_region = 5,
                                 mc.cores = 18
  )
  print(Sys.time() - time_start)
  
  ####
  
  regions_all_chr <- dap_reduced[[1]]$regions
  coverage_all_chr <- dap_reduced[[1]]$coverage
  methylation_all_chr <- dap_reduced[[1]]$methylation
  for (i in 2:length(dap_reduced)) {
    regions_all_chr <- rbind(regions_all_chr, dap_reduced[[i]]$regions)
    coverage_all_chr <- rbind(coverage_all_chr, dap_reduced[[i]]$coverage)
    methylation_all_chr <- rbind(methylation_all_chr, dap_reduced[[i]]$methylation) 
  }
  
  dim(coverage_all_chr)
  
  #save out all this work
  write_rds(regions_all_chr, file = paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_regions_oi.rds"))
  write_rds(coverage_all_chr, file = paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_coverage_regions_oi.rds"))
  write_rds(methylation_all_chr, file = paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_methylation_regions_oi.rds"))
}
