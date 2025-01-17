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
  n_cores = 4
  chr="1"
  batch = 1
} else {
  args <- commandArgs(trailingOnly=TRUE)
  chr <- args[1]
  if(chr==39){chr="X"}
  batch <- args[2]
  n_cores = 12
}

if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
if(dir.exists("pqlseq_results") == F){dir.create("pqlseq_results")}

library_list <- c(
  "PQLseq2",
  "VariantAnnotation"
)
lapply(library_list, require, character.only = TRUE)

metaData <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/metaData_samples/all_precision_metaData.rds")
metaData <- metaData[metaData$use_for_baseline_precision_clock == "yes",]
metaData<-metaData[order(metaData$dog_id),]

library(genio)
print(paste("chr", chr))

plink=FALSE
if(plink==TRUE){
  plink_data <- read_plink(paste0("/scratch/bmarine2/meQTL-241101-GEMMA/gemma_out/filtered_genotypes_chr",chr,".bed"))
  plink_in <- plink_data$X
  snp_bims <- data.frame(plink_data$bim)
  snp_bims$chr_snp <- paste0(chr, "_", snp_bims$pos)
  rownames(plink_in) <- snp_bims$chr_snp
  plink_in <- t(plink_in)
  plink_in <- (plink_in[,colSums(plink_in == 2) < .9*nrow(plink_in) & colSums(plink_in == 1) < .9*nrow(plink_in) & colSums(plink_in == 0) < .9*nrow(plink_in)])

  colnames(plink_in) <- gsub(paste0(chr,"_"), "", colnames(plink_in))

  if(chr == "X"){
    colnames(plink_in) <- gsub(paste0("39_"), "", colnames(plink_in))
  }

  rownames(plink_in) <- gsub(".canfam4", "", rownames(plink_in))
  rm(plink_data)

  # Combine the results into a dataframe
  plink_in<- data.frame(t(unique(t(plink_in))))
  colnames(plink_in) <- gsub("X", "", colnames(plink_in))
} else {
  DS_scores <- read.table(paste0("/scratch/bmarine2/vcf_fun-250108/matrix_outputs/combined_chr",chr,"_matrix.txt"))
  rownames(DS_scores) <- DS_scores[,1]
  colnames(DS_scores) <- DS_scores[1,]
  DS_scores <- DS_scores[2:nrow(DS_scores), 2:ncol(DS_scores)]

  rownames(DS_scores) <- gsub(":","_",rownames(DS_scores))

  INFO <- read.table(paste0("/scratch/bmarine2/vcf_fun-250108/info_outputs/combined_chr",chr,"_info_values.txt"))
  colnames(INFO) <- c("chr_snp", "INFO")
  INFO <- INFO[INFO$INFO > 0.5,]
  DS_scores <- DS_scores[rownames(DS_scores) %in% INFO$chr_snp,]

  filter_variab=TRUE
  if(filter_variab==T){
    rownames_og <- rownames(DS_scores)
    DS_scores <- as.data.frame(lapply(DS_scores, type.convert, as.is = TRUE))
    colnames(DS_scores) <- gsub("X","",colnames(DS_scores))
    rownames(DS_scores) <- rownames_og
    row_mean <- rowMeans(DS_scores, na.rm = TRUE)

    # Set a threshold for adequate variability (e.g., 0.01)
    thresholdlow <- 0.2
    thresholdhigh <- 2
    # Filter rows with SD above the threshold
    filtered_df <- DS_scores[row_mean > thresholdlow & row_mean < thresholdhigh, ]
  }

  plink_in <- t(filtered_df)
  rm(DS_scores)
}


groups <- split(1:ncol(plink_in),
                cut(1:ncol(plink_in), breaks = 1000, labels = FALSE)
)

plink_in <- plink_in[,groups[[as.numeric(batch)]]]
colnames(plink_in) <- gsub(".*_","",colnames(plink_in))

plink_in_ps<-colnames(plink_in)

RelatednessMatrix_ <- readRDS(file = "/scratch/bmarine2/BASELINE-PRECISION-240820/metaData_samples/240821-baseline_GRM.rds")

regions <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_regions_1n_coverage.rds")
cov =  readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_coverage_regions_1n_coverage.rds")
meth =  readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_methylation_regions_1n_coverage.rds")
regions_oi <- regions[regions$chr == paste0("chr", chr),]
cov_oi <- cov[startsWith(prefix = paste0("chr", chr, "_"), x = rownames(cov)),]
meth_oi <- meth[startsWith(prefix = paste0("chr", chr, "_"), x = rownames(meth)),]
ps_func <- function(ps, regions_oi, cov_oi, meth_oi, distance=50000){
  if(exists("cov_append")){rm(cov_append)}
  if(exists("meth_append")){rm(meth_append)}

  if(file.exists(paste0("coverage_meth_snp_paired/snp_cov_chr", chr, "_", ps, ".tsv"))){return(NULL)}
  # print(ps)
  ps<-as.numeric(ps)
  rangesB <- IRanges::IRanges(regions_oi$start-distance, regions_oi$end+distance)
  ov <- GenomicRanges::countOverlaps(rangesB, IRanges::IRanges(ps), type="any")>0
  ## this ps list was pre-filtered... so this shouldnt happen if(all(ov==F)){return(NA)}

  hit_ <- regions_oi[ov,]
  hit_$snp <- rep(ps, nrow(hit_))
  if(nrow(hit_) == 0){return(NULL)}
  hit_$new_rowname <- paste0(rownames(hit_), "_", hit_$snp)
  for (rowname in rownames(hit_)){

    meth_append_ <- meth_oi[rownames(meth_oi) == rowname,,drop=F]
    cov_append_ <- cov_oi[rownames(cov_oi) == rowname,,drop=F]

    if(exists("cov_append") == F) {
      rownames(meth_append_) <- rownames(cov_append_) <- hit_$new_rowname[rownames(hit_) == rowname]
      cov_append <- cov_append_
      meth_append <- meth_append_
    }else{
      rownames_og <- rownames(cov_append)
      cov_append <- rbind(cov_append, cov_append_)
      meth_append <- rbind(meth_append, meth_append_)
      rownames(meth_append) <- rownames(cov_append) <- c(rownames_og, hit_$new_rowname[rownames(hit_) == rowname])
    }
  }
  list_res <- list(cov_append, meth_append)
  names(list_res) <- c("cov", "meth")

  write.table(list_res[["cov"]], paste0("coverage_meth_snp_paired/snp_cov_chr", chr, "_", ps, ".tsv"))
  write.table(list_res[['meth']], paste0("coverage_meth_snp_paired/snp_meth_chr", chr, "_", ps, ".tsv"))

  return(list_res)
}

print("pairing snps with regions")
res_list <- parallel::mclapply(unique(plink_in_ps),
                               ps_func,
                               regions_oi,cov_oi, meth_oi,
                               mc.cores=n_cores
)
#gonna read in from file for the next thing, so to save mem:
rm(res_list)
rm(meth)
rm(cov)

pqlseq_snp_fnct <- function (ps, chr, plink_in, metaData, RelatednessMatrix_){

  if(! file.exists(paste0("pqlseq_results/pqlseq_res", chr, "_", ps, ".tsv"))){

    if(! file.exists(paste0("coverage_meth_snp_paired/snp_cov_chr", chr, "_", ps, ".tsv"))){return(NULL)}

    cov_oi <- read.table(paste0("coverage_meth_snp_paired/snp_cov_chr", chr, "_", ps, ".tsv"))
    meth_oi <- read.table(paste0("coverage_meth_snp_paired/snp_meth_chr", chr, "_", ps, ".tsv"))

    tmp <- as.data.frame(do.call("rbind", stringr::str_split(rownames(cov_oi), "_")))
    colnames(tmp) <- c("chr", "start", "end", "snp")

    metaData <- metaData[metaData$dog_id %in% rownames(plink_in),]
    metaData <- metaData[order(metaData$dog_id),]
    colnames_og <- colnames(metaData)
    plink_in <- plink_in[order(rownames(plink_in)),]

    if(!(all(metaData$dog_id == rownames(plink_in)))){stop("sort your shtuff")}

    metaData[,paste(ps)] <- (plink_in[,c(paste(ps)),drop=F])
    metaData$SNP <- metaData[,paste(ps)]
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
    # print(pqlseq_res_snp_)
    return(print(paste0("done ", chr, "_", ps)))
  }
}

print("pqlseq2")
pqlseq_list_snps <- lapply(plink_in_ps,
                           pqlseq_snp_fnct,
                           chr, plink_in, metaData,RelatednessMatrix_
)
print((Sys.time() - time_start))


