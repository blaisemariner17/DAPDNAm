#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript
#### Blaise Mariner 
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-01-12
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project
rm(list = ls())
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  n_cores = 1
} else {n_cores = 20}

library_list <- c(
  "ggplot2",
  "svglite",
  "ggExtra",
  "ggtext",
  "GenomicRanges",
  "GenomicFeatures",
  "tidyverse",
  "corrplot",
  "glmnet",
  "PQLseq",
  "IRanges",
  "patchwork"
)
lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 12, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      
                      axis.text = element_markdown(family = "sans", size = 14, color="black"),axis.title.y = element_markdown(family = "sans", size = 20),   
                      axis.title.x = element_markdown(family = "sans", size = 20),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

#1.5 hrs i
metaData <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/metaData_samples/all_precision_metaData.rds")
metaData <- metaData[metaData$use_for_baseline_precision_clock == "yes",]
metaData<-metaData[order(metaData$dog_id),]

# files_oi <- list.files("./coverage_meth_snp_paired/",  pattern = "snp_cov_chr1_")
# region_snp_all <- 0
# i=1
# for (file in files_oi){
#   region_snp_all <- region_snp_all + length(rownames(read.table(paste0("./coverage_meth_snp_paired/",file))))
#   print(region_snp_all)
#   print(i/length(files_oi))
# }


load("in_progress.RData")

# install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
## https://marinalearning.netlify.app/2021/03/22/setting-up-multivariable-mendelian-randomization-analysis/
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)
# Install and load required package
if (!requireNamespace("MendelianRandomization", quietly = TRUE)) {
  install.packages("MendelianRandomization")
}
library(MendelianRandomization)
library(parameters)

chr = 1

for (chr in 1:38){
  print(chr)
  
  pqlseq_res_oi <- pqlseq_res[pqlseq_res$chr == paste0("chr", chr),]
  
  plink_in <- MultiPhen::read.plink(paste0("/scratch/bmarine2/meQTL-241101-GEMMA/gemma_out/filtered_genotypes_chr", chr))
  plink_in<-plink_in[order(rownames(plink_in)),]
  metaData <- metaData[metaData$dog_id %in% rownames(plink_in),]
  metaData <- metaData[order(metaData$lid_pid),]
  
  cov <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_coverage_regions_1n_coverage.rds")
  meth <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_methylation_regions_1n_coverage.rds")
  
  cov <- cov[startsWith(rownames(cov), paste0("chr",chr,"_")),]
  meth <- meth[startsWith(rownames(meth), paste0("chr",chr,"_")),]
  
  cov <- cov[,colnames(cov) %in% metaData$lid_pid]
  meth <- meth[,colnames(meth) %in% metaData$lid_pid]
  colnames(cov) <- colnames(meth) <- metaData$dog_id
  
  pqlseq_res_oi <- pqlseq_res_oi[order(pqlseq_res_oi$pvalueAge_at_sample),]
  if(nrow(pqlseq_res_oi) == 0){next}
  
  snp <- pqlseq_res_oi$snp[1]
  region_snp <- rownames(pqlseq_res_oi)[1]
  region <- pqlseq_res_oi$region[1]
  
  cov_oi <- cov[region,order(colnames(cov)), drop = F]
  meth_oi <- meth[region,order(colnames(meth)), drop = F]
  perc_meth_oi <- t(meth_oi / cov_oi)
  perc_meth_oi <- perc_meth_oi[,order(colnames(perc_meth_oi))]
  
  # perc_meth <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/perc_meth_imputed/methyLImp2_perc_meth_1n_imputed-241017.rds")
  # perc_meth <- perc_meth[,colnames(perc_meth) %in% metaData$lid_pid]
  # perc_meth <- perc_meth[,order(colnames(perc_meth))]
  # colnames(perc_meth) <- metaData$dog_id
  # perc_meth_oi<-perc_meth[region,order(colnames(perc_meth))]
  
  plink_oi <- (plink_in[,c(paste0(chr, "_", snp)),drop=F])
  a_freq <- ((table(plink_oi))[1] + (table(plink_oi)[2]/2)) / sum((table(plink_oi)))
  metaData <- metaData[order(metaData$dog_id),]
  if (!all(metaData$dog_id == rownames(plink_oi))){stop("fix plink and metadata dog_ids")}
  if (!all(metaData$dog_id == colnames(perc_meth_oi))){stop("fix plink and metadata dog_ids")}
  
  metaData$SNP <- as.data.frame(plink_oi)[,1]
  metaData$region_pm <- perc_meth_oi
  
  beta_meth_x <- pqlseq_res_oi$betaSNP    # SNP effect on methylation from the binomial mixed model
  se_meth_x <- pqlseq_res_oi$se_betaSNP         # Standard error for each SNP effect

  beta_meth_y <- pqlseq_res_oi$betaAge_at_sample    # SNP effect on methylation from the binomial mixed model
  se_meth_y <- pqlseq_res_oi$se_betaAge_at_sample         # Standard error for each SNP effect
  
  
  # Step 2: Set up the data for MR
  # Assuming we're running MR for each SNP individually
  # This example illustrates the IVW method for each SNP effect estimate
  
  # Loop through each SNP in `pqlseq_res_oi` and perform MR
  mr_results <- lapply(1:nrow(pqlseq_res_oi), function(i) {
    # Create MR input for each SNP
    mr_data <- mr_input(bx = beta_meth_x[i], bxse = se_meth_x[i], by = beta_meth_y[i], byse = se_meth_y[i])
    
    # Perform Mendelian Randomization with IVW method
    mr_ivw_result <- mr_ivw(mr_data)
    
    # Store results in a list
    list(
      SNP = rownames(pqlseq_res_oi)[i],
      MR_estimate = mr_ivw_result$Estimate,
      MR_se = mr_ivw_result$StdError,
      MR_pval = mr_ivw_result$Pvalue
    )
  })
  
  # Convert results to a data frame for easier viewing
  if(chr == 1){mr_results_df <- do.call(rbind, lapply(mr_results, as.data.frame))} else{
    mr_results_df <- rbind(mr_results_df, do.call(rbind, lapply(mr_results, as.data.frame)))
  }
}

hist(mr_results_df$MR_estimate)
hist(mr_results_df$MR_pval)

mr_results_df$MR_padj <- p.adjust(mr_results_df$MR_pval, method = 'fdr')
mr_results_df <- mr_results_df[order(mr_results_df$MR_padj),]

mr_results_df$chr_snp <- gsub("\\..*", "", mr_results_df$SNP)
mr_results_df$region <- gsub(".*\\.", "", mr_results_df$SNP)
mr_results_df$region <- gsub("^(([^_]*_){3}).*", "\\1", mr_results_df$region)
mr_results_df$region <- sub("_$", "", mr_results_df$region)

head(mr_results_df)

region_metaData <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/metadata_regions/metaData_regions_maxgap_1n_250.rds")
col_oi = c(
  "Promoter", #"gene_bool",
  "exon", "intron",
  "CpG_island", "CpG_shore", "CpG_shelf",
  # "ChrSt_promoter",
  "ChrSt_quies","ChrSt_heterochromatin",
  # "ChrSt_polycomb",
  "ChrSt_enhancer",
  # "DNA.transposon", "Retrotransposon",
  "TE"
)

fit_fdr_snp <- mr_results_df[mr_results_df$MR_padj < 0.1,]

region_metaData$fdr_snp <- "not significant"
region_metaData$fdr_snp[region_metaData$region %in% fit_fdr_snp$region] <- "SNP influence"
region_metaData <- region_metaData[region_metaData$region %in% mr_results_df$region,]

i=1
for (col in c(col_oi)){
  a <- nrow(region_metaData[region_metaData[,paste(col)] == 1  & region_metaData$fdr_snp == "SNP influence",])
  b <- nrow(region_metaData[region_metaData[,paste(col)] == 1 & region_metaData$fdr_snp == "not significant",])
  c <- nrow(region_metaData[(!region_metaData[,paste(col)] == 1)  & region_metaData$fdr_snp == "SNP influence",])
  d <- nrow(region_metaData[(!region_metaData[,paste(col)] == 1) & region_metaData$fdr_snp == "not significant",])
  
  data_oi <- matrix(c(a,b,c,d),
                    nrow=2, ncol=2, byrow=TRUE)
  
  odds_ratio <- epitools::oddsratio(data_oi)
  res_ <- data.frame("class" = col, "pval" = round(odds_ratio$p.value[2,2], digits = 10), "odds_ratio_log2" = round(log2(odds_ratio$measure[2,1]), digits = 10),
                     "lower" = round(log2(odds_ratio$measure[2,2]), digits = 10), "upper" = round(log2(odds_ratio$measure[2,3]), digits = 10))
  if (i == 1){
    res <- res_
    i = 2
  } else {
    res <- rbind(res, res_)
  }
}
res$padj <- p.adjust(res$pval, method = "BH")
res <- res[order(res$pval),]
res


