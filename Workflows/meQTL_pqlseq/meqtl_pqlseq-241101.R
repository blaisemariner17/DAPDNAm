#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript
#### Blaise Mariner 
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-01-12
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project
rm(list = ls())  # Clear the workspace
time_start <- Sys.time()

# Determine if running in RStudio; set the number of cores accordingly
if (isRStudio <- Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set working directory
  n_cores = 4  # Single-core for RStudio
  chr=32
  batch = 55
} else {
  args <- commandArgs(trailingOnly=TRUE)
  chr <- args[1]
  batch <- args[2]
  n_cores = 20
}

if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
if(dir.exists("pqlseq_results") == F){dir.create("pqlseq_results")}

library_list <- c(
  "PQLseq",
  "VariantAnnotation"
)
lapply(library_list, require, character.only = TRUE)

metaData <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/metaData_samples/all_precision_metaData.rds")
metaData <- metaData[metaData$use_for_baseline_precision_clock == "yes",]
metaData<-metaData[order(metaData$dog_id),]


print(paste("chr", chr))
# plink_in <- MultiPhen::read.plink("/scratch/vsohrab/dap_genetic_set_2023/DogAgingProject_2023_N-7627_canfam4_gp-0.70_biallelic")
plink_in <- MultiPhen::read.plink(paste0("/scratch/bmarine2/meQTL-241101-GEMMA/gemma_out/filtered_genotypes_chr", chr))
#filter for variable snps in your samples-- plink should have already done this
plink_in <- (plink_in[,colSums(plink_in == 2) < .9*nrow(plink_in) & colSums(plink_in == 1) < .9*nrow(plink_in)])
colnames(plink_in) <- gsub(paste0(chr,"_"), "", colnames(plink_in))

remove_duplicates <- function(col_name, data) {
  grep_ <- grep(col_name, colnames(data))
  if (length(grep_) > 1) {
    # Return only the first column
    return(data[, grep_[1], drop = FALSE])
  } else {
    # Return the single column if no duplicates
    return(data[, grep_, drop = FALSE])
  }
}

# Parallel processing
result_list <- parallel::mclapply(unique(colnames(plink_in)), function(col) {
  remove_duplicates(col, plink_in)
}, mc.cores = n_cores) # Use all cores but one

# Combine the results into a dataframe
plink_in<- do.call(cbind, result_list)

i_list <- 1:ncol(plink_in)
i_for_par <- split(i_list, ceiling(seq_along(i_list)/100))
plink_in <- plink_in[,i_for_par[[as.numeric(batch)]]]

plink_in_ps<-colnames(plink_in)

RelatednessMatrix_ <- readRDS(file = "/scratch/bmarine2/BASELINE-PRECISION-240820/metaData_samples/240821-baseline_GRM.rds")

# regions <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_regions_1n_coverage.rds")
# cov =  readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_coverage_regions_1n_coverage.rds")
# meth =  readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_methylation_regions_1n_coverage.rds")
# regions_oi <- regions[regions$chr == paste0("chr", chr),]
# cov_oi <- cov[startsWith(prefix = paste0("chr", chr, "_"), x = rownames(cov)),]
# meth_oi <- meth[startsWith(prefix = paste0("chr", chr, "_"), x = rownames(meth)),]
# ps_func <- function(ps, regions_oi, cov_oi, meth_oi, distance=50000){
#   if(exists("cov_append")){rm(cov_append)}
#   if(exists("meth_append")){rm(meth_append)}
#   # print(ps)
#   ps<-as.numeric(ps)
#   rangesB <- IRanges::IRanges(regions_oi$start-distance, regions_oi$end+distance)
#   ov <- GenomicRanges::countOverlaps(rangesB, IRanges::IRanges(ps), type="any")>0
#  ## this ps list was pre-filtered... so this shouldnt happen if(all(ov==F)){return(NA)}
#   hit_ <- regions_oi[ov,]
#   hit_$snp <- rep(ps, nrow(hit_))
#   hit_$new_rowname <- paste0(rownames(hit_), "_", hit_$snp)
#   for (rowname in rownames(hit_)){
# 
#     meth_append_ <- meth_oi[rownames(meth_oi) == rowname,,drop=F]
#     cov_append_ <- cov_oi[rownames(cov_oi) == rowname,,drop=F]
# 
#     if(exists("cov_append") == F) {
#       rownames(meth_append_) <- rownames(cov_append_) <- hit_$new_rowname[rownames(hit_) == rowname]
#       cov_append <- cov_append_
#       meth_append <- meth_append_
#     }else{
#       rownames_og <- rownames(cov_append)
#       cov_append <- rbind(cov_append, cov_append_)
#       meth_append <- rbind(meth_append, meth_append_)
#       rownames(meth_append) <- rownames(cov_append) <- c(rownames_og, hit_$new_rowname[rownames(hit_) == rowname])
#     }
#   }
#   list_res <- list(cov_append, meth_append)
#   names(list_res) <- c("cov", "meth")
# 
#   write.table(list_res[["cov"]], paste0("coverage_meth_snp_paired/snp_cov_chr", chr, "_", ps, ".tsv"))
#   write.table(list_res[['meth']], paste0("coverage_meth_snp_paired/snp_meth_chr", chr, "_", ps, ".tsv"))
# 
#   return(list_res)
# }
# res_list <- parallel::mclapply(unique(plink_in_ps),
#                                ps_func,
#                                regions_oi,cov_oi, meth_oi,
#                                mc.cores=n_cores
# )
# #gonna read in from file for the next thing, so to save mem:
# rm(res_list)

pqlseq_snp_fnct <- function (ps, chr, plink_in, metaData, RelatednessMatrix_){
  
  if(! file.exists(paste0("pqlseq_results/pqlseq_res", chr, "_", ps, ".tsv"))){
    cov_oi <- read.table(paste0("coverage_meth_snp_paired/snp_cov_chr", chr, "_", ps, ".tsv"))
    meth_oi <- read.table(paste0("coverage_meth_snp_paired/snp_meth_chr", chr, "_", ps, ".tsv"))
    
    tmp <- as.data.frame(do.call("rbind", stringr::str_split(rownames(cov_oi), "_")))
    colnames(tmp) <- c("chr", "start", "end", "snp")
    
    metaData <- metaData[metaData$dog_id %in% rownames(plink_in),]
    metaData <- metaData[order(metaData$dog_id),]
    plink_in <- plink_in[order(rownames(plink_in)),]
    
    if(!(all(metaData$dog_id == rownames(plink_in)))){stop("sort your shtuff")}
    
    metaData$SNP <- (plink_in[,c(paste(ps)),drop=F])
    
    cov_oi <- cov_oi[,order(colnames(cov_oi))]
    meth_oi <- meth_oi[,order(colnames(meth_oi))]
    metaData <- metaData[order(metaData$lid_pid),]
    
    design <- model.matrix(
      ~SNP + Age_at_sample + Sex_bool, data = metaData
    )
    
    fit_colname <- colnames(design)[2]
    new_colname <- paste0(fit_colname)
    
    pheno <- (design[,paste(fit_colname),drop=F])
    covariates <- as.matrix(design[,colnames(design)[colnames(design) != fit_colname & colnames(design) != "(Intercept)"]])
    
    meth_oi <- meth_oi[,colnames(meth_oi) %in% metaData$lid_pid, drop=F]
    cov_oi <- cov_oi[,colnames(cov_oi) %in% metaData$lid_pid, drop=F]
    
    RelatednessMatrix_<- RelatednessMatrix_[rownames(RelatednessMatrix_) %in% colnames(cov_oi), colnames(RelatednessMatrix_) %in% colnames(cov_oi)]
    
    pqlseq_res_snp_ = PQLseq2::pqlseq2(Y=meth_oi,
                                       x=pheno,
                                       K=RelatednessMatrix_,
                                       lib_size=cov_oi,
                                       model="BMM",
                                       W = covariates,
                                      ncores = n_cores)
    
    # pqlseq_res_snp_ = PQLseq::pqlseq(RawCountDataSet=meth_oi[1:2,],
    #                                     Phenotypes=pheno,
    #                                     RelatednessMatrix=RelatednessMatrix_,
    #                                     LibSize=cov_oi[1:2,],
    #                                     fit.model="BMM",
    #                                     Covariates = covariates,
    #                                     numCore = n_cores)
    
    pqlseq_res_snp_ <- pqlseq_res_snp_[order(rownames(pqlseq_res_snp_)),]
    colnames(pqlseq_res_snp_) <- paste0(colnames(pqlseq_res_snp_), new_colname)

    write.table(pqlseq_res_snp_, paste0("pqlseq_results/pqlseq_res", chr, "_", ps, ".tsv"))
    print(pqlseq_res_snp_)
    return(print(paste0("done ", chr, "_", ps)))
  }
}

pqlseq_list_snps <- lapply(plink_in_ps,
                                      pqlseq_snp_fnct,
                                      chr, plink_in, metaData,RelatednessMatrix_
                  )
print((Sys.time() - time_start))


