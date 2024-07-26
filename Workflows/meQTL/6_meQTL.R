#!/usr/bin/Rscript
## this script is supposed to do the same thing as : https://github.com/Cec701/genetic_architecture_DNAm_rhesus/blob/main/IMAGE_meQTL/chrn_573_merge_filter.R
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

## example output of script:
# load("/scratch/ccosta8/RRBS/ASM_all/merge_filter_chr18_573.RData")

#Clock manuscript figure generation - 05/14/24

#load libraries needed
library(tidyverse)
library(ggplot2)
library(ggtext)
library(ggExtra)
library(patchwork)
library(IMAGE)

compile_asms=F
if (compile_asms){
  i=1
  for(file_oi in list.files("ASM/")){
    print(file_oi)
    asm_oi <- read.table(paste0("ASM/", file_oi), header = T)
    asm_oi$lid_pid <- gsub(".CGmap.PASS.DP5.ASM","",file_oi)
    if(i==1){asm<-asm_oi;i=2}else{asm<-rbind(asm,asm_oi)}
  }
  write.table(asm,"ALL_ASM.ASM")
}


chrs <- paste0("chr", c(1:38))

getSNPs <- function(chr) {
  print(chr)
  NUMIDS <- length(list.files("ASM/"))

  asm_all <- data.frame(read.table("ALL_ASM.ASM", header = T))

  asm_chr <- asm_all[asm_all$Chr == chr,]
  rm(asm_all)

  #remove the cpgs that are in less than 10% of individuals (Fan et al., 2019)
  asm_chr <- asm_chr[asm_chr$C_Pos %in% names(table(asm_chr$C_Pos)[table(asm_chr$C_Pos) > .1*NUMIDS]),]

  #remove the snp cpg instances found with >5 total reads (Fan et al., 2019)
  asm_chr$Allele1_reads <- as.numeric(gsub("-.*","", asm_chr$Allele1_linked_C)) + as.numeric(gsub(".*-","", asm_chr$Allele1_linked_C))
  asm_chr$Allele2_reads <- as.numeric(gsub("-.*","", asm_chr$Allele2_linked_C)) + as.numeric(gsub(".*-","", asm_chr$Allele2_linked_C))
  asm_chr$total_reads <- asm_chr$Allele1_reads + asm_chr$Allele2_reads
  asm_chr <- asm_chr[asm_chr$total_reads >=5,]

  #remove the perpetually methylated snp cpgs >90% or <10% (Fan et al., 2019)
  asm_chr$Allele1_meth <- as.numeric(gsub("-.*","", asm_chr$Allele1_linked_C))
  asm_chr$Allele2_meth <- as.numeric(gsub("-.*","", asm_chr$Allele2_linked_C))

  asm_chr$Allele1_pmeth <- asm_chr$Allele1_meth / asm_chr$Allele1_reads
  asm_chr$Allele2_pmeth <- asm_chr$Allele2_meth / asm_chr$Allele2_reads

  asm_chr <- asm_chr[asm_chr$Allele1_pmeth < 0.9 | asm_chr$Allele1_pmeth > 0.1 | asm_chr$Allele2_pmeth < 0.9 | asm_chr$Allele2 > 0.1,]

 # for (snp in unique(asm_chr$SNP_Pos)){
 #   mean_pmeth1 <- mean(asm_chr$Allele1_pmeth[asm_chr$SNP_Pos==snp])
 #   mean_pmeth2 <- mean(asm_chr$Allele2_pmeth[asm_chr$SNP_Pos==snp])
 #   if((mean_pmeth1 > 0.90 & mean_pmeth2 > 0.90 ) | (mean_pmeth1 < 0.1 & mean_pmeth2 < 0.1)){asm_chr<-asm_chr[asm_chr$SNP_Pos != snp,]}
 # }

  # require at least 5% mean allele frequency up to here (Fan et al., 2019)
  asm_chr <- asm_chr[asm_chr$SNP_Pos %in% names(table(asm_chr$SNP_Pos)[table(asm_chr$SNP_Pos) > .05*NUMIDS]),]

  #label the geno of the snp
  # and their "geno" variable is also infuriating but 0 means the reference allele matches both allele 1 and allele 2, 2 means neither match, and 1 means else....
  asm_chr$geno <- 1
  asm_chr$geno[asm_chr$Ref == asm_chr$Allele1 & asm_chr$Allele1 == asm_chr$Allele2] <- 0
  asm_chr$geno[asm_chr$Ref != asm_chr$Allele1 & asm_chr$Allele1 != asm_chr$Allele2] <- 2

  return(asm_chr)
}

getSNPs_list <- parallel::mclapply(chrs,
                                   getSNPs,
                                   mc.cores = 20
                                   )

save.image(file = "meQTL_filtering/merge_filter.RData")


##########################      meqtl       ################################################3
# we have to prepare the data because the above function gets everything done faster but now the data structure has to be corrected
# note their non-descriptive variables: r = read count, y = methylated count. 1 = allele 1 and 2 = allele 2. so y1 for example is the methylated count for allele 1

format_for_image <- function(getSNPs){
  rownames(getSNPs) <- paste(getSNPs$Chr, getSNPs$SNP_Pos, getSNPs$C_Pos, getSNPs$lid_pid, sep='_')
  getSNPs <- getSNPs[order(rownames(getSNPs)),]
  getSNPs <- getSNPs[order(getSNPs$lid_pid),]

  lid_pids <- unique(getSNPs$lid_pid)
  SNP_CPG <- unique(gsub("_L.*", "", rownames(getSNPs)))

  r <- data.frame(matrix(0,ncol=length(lid_pids),nrow=length(SNP_CPG)))
  colnames(r) <- lid_pids
  rownames(r) <- SNP_CPG
  geno <- y2 <- y1 <- y <- r2 <- r1 <- r

  for (snp_cpg in SNP_CPG){
    snpcpg_oi <- getSNPs[grepl(snp_cpg, rownames(getSNPs)),]

    lid_pids_with_snp <- unique(snpcpg_oi$lid_pid)
    lid_pids_not_with_snp <- lid_pids[!lid_pids %in% lid_pids_with_snp]


    appending_lid_pids <- data.frame(matrix(0, ncol = ncol(snpcpg_oi), nrow = length(lid_pids_not_with_snp)))
    colnames(appending_lid_pids) = colnames(snpcpg_oi)
    appending_lid_pids$lid_pid <- lid_pids_not_with_snp

    snpcpg_oi <- rbind(snpcpg_oi, appending_lid_pids)

    snpcpg_oi <- snpcpg_oi[order(rownames(snpcpg_oi)),]
    snpcpg_oi <- snpcpg_oi[order(snpcpg_oi$lid_pid),]

    genotype <- snpcpg_oi$geno
    total_reads <- snpcpg_oi$total_reads
    allele1_reads <- snpcpg_oi$Allele1_reads
    allele2_reads <- snpcpg_oi$Allele2_reads
    total_meth_reads <- snpcpg_oi$Allele1_meth + snpcpg_oi$Allele2_meth
    allele1_meth_reads <- snpcpg_oi$Allele1_meth
    allele2_meth_reads <- snpcpg_oi$Allele2_meth

    geno[rownames(geno) == snp_cpg,] <- genotype
    r[rownames(r) == snp_cpg,] <- total_reads
    r1[rownames(r1) == snp_cpg,] <- allele1_reads
    r2[rownames(r2) == snp_cpg,] <- allele2_reads
    y[rownames(y) == snp_cpg,] <- total_meth_reads
    y1[rownames(y1) == snp_cpg,] <- allele1_meth_reads
    y2[rownames(y2) == snp_cpg,] <- allele2_meth_reads
  }

  data_res <- list()
  data_res[['geno']] <- geno
  data_res[['r']] <- r
  data_res[['r1']] <- r1
  data_res[['r2']] <- r2
  data_res[['y']] <- y
  data_res[['y1']] <- y1
  data_res[['y2']] <- y2
  return(data_res)
}

IMAGE_format_list <- parallel::mclapply(getSNPs_list,
                                        format_for_image,
                                   mc.cores = 20
)

save.image(file = "meQTL_filtering/merge_filter.RData")

IMAGE_format <- do.call("rbind", IMAGE_format_list)

IMAGE_format_data <- IMAGE_format[names(IMAGE_format) %in% c("r","y","r1", "r2","y1","y2")]
IMAGE_format_geno<-IMAGE_format[["geno"]]

RelatednessMatrix <- readRDS(file = "../metadata_samples/240415-p123_GRM.rds")

IMAGE::image(geno = IMAGE_format_geno,
             data = IMAGE_format_data,
             K = RelatednessMatrix,
             numCore=24)

