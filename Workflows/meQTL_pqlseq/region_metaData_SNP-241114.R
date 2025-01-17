#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript
#### Blaise Mariner
### this script generates the region metaData such that there is a way to better understand the reigons that are differentially methylated with different demographics
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-02-13
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project

#### GLOBALS ####
# this scrip requires at least 12 cores
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
  "Biostrings",
  "patchwork",
  "stringr"
)

lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 24, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      axis.text = element_text(family = "sans", size = 18, color="black"),       axis.title = element_text(family = "sans", size = 20, color="black"),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  n_cores = 4
  batch = 1000
} else {
  args <- commandArgs(trailingOnly=TRUE)

  batch <- args[1]
  n_cores = 12
}
print(getwd())

gtf <- rtracklayer::import('/home/bmarine2/GENOME-ANNOTATION-FILE/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf')
gtf_cpgisl <- read.csv("/home/bmarine2/GENOME-ANNOTATION-FILE/cpgIslandExt.txt", header = F, sep = "\t")
bed_file <- "/home/bmarine2/GENOME-ANNOTATION-FILE/SP_13_dense.bed"
chromatin_states_annotation = rtracklayer::import(bed_file, format = "bed")
chromatin_states_annotation$name <- as.character(chromatin_states_annotation$name)
chromatin_states_annotation$name[chromatin_states_annotation$name == '1' | chromatin_states_annotation$name == '2'|
                                   chromatin_states_annotation$name == '3' | chromatin_states_annotation$name == '4'] <- "ChrSt_promoter"
chromatin_states_annotation$name[chromatin_states_annotation$name == '5' | chromatin_states_annotation$name == '6'|
                                   chromatin_states_annotation$name == '7'] <- "ChrSt_enhancer"
chromatin_states_annotation$name[chromatin_states_annotation$name == '8' | chromatin_states_annotation$name == "9"] <- "ChrSt_polycomb"
chromatin_states_annotation$name[chromatin_states_annotation$name == '10' | chromatin_states_annotation$name == "11"|
                                   chromatin_states_annotation$name == "12"] <- "ChrSt_heterochromatin"
chromatin_states_annotation$name[chromatin_states_annotation$name == "13"] <- "ChrSt_quies"

gtf_trans <- read.csv("/home/bmarine2/GENOME-ANNOTATION-FILE/rmsk.txt", header = F, sep = "\t")
colnames(gtf_trans) <- c('bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns'	,
                         'seqnames', 'start', 'end', 'genoLeft', 'strand'	,
                         'repName', 'class', 'repFamily', 'repStart', 'repEnd',	'repLeft',	'id')
gtf_trans <- gtf_trans[gtf_trans$seqnames %in% paste0("chr", c(1:38, "X")),]

# gtf_trans$length<-gtf_trans$end-gtf_trans$start
# for_gg <- gtf_trans[gtf_trans$class == "LINE",]
# for_gg <- for_gg[grepl("L1M", for_gg$repName),]
# ggplot(for_gg,aes(x=repName, y = length))+
#   geom_boxplot()+theme_blaise+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


transposon_annotation <- gtf_trans[gtf_trans$class %in% c("LINE", "SINE", "LTR", "DNA",
                                                          "Simple_repeat","Low_complexity", "rRNA", "tRNA",
                                                          "Satellite","Unknown"
),]

dim(transposon_annotation)

rm(gtf_trans)

# .rs.restartR()
# remove.packages("DAPDNAm")
if (! "DAPDNAm" %in% installed.packages()){devtools::install_github("blaisemariner17/DAPDNAm", force = TRUE)}
library("DAPDNAm")

if(dir.exists("metaData_SNPs")==F){dir.create("metaData_SNPs")}

# regions <- read_rds(paste0("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_regions_1n_coverage.rds"))

if(TRUE){
  library(genio)
  for (chr in c(1:38,"X")){
    if(file.exists(paste0("metaData_SNPs/chr",chr,"_SNPs.tsv"))){next}

    print(paste0("chr", chr))

    plink_data <- read_plink(paste0("/scratch/bmarine2/meQTL-241101-GEMMA/gemma_out/filtered_genotypes_chr",chr,".bed"))
    plink_in <- plink_data$X
    snp_bims <- unique(data.frame(plink_data$bim))

    # snp_bims$chr_snp <- paste0(chr, "_", snp_bims$pos)
    # rownames(plink_in) <- snp_bims$chr_snp
    # plink_in <- t(plink_in)
    #
    # colnames(plink_in) <- gsub(paste0(chr,"_"), "", colnames(plink_in))
    #
    # if(chr == "X"){
    #   colnames(plink_in) <- gsub(paste0("39_"), "", colnames(plink_in))
    # }
    #
    # rownames(plink_in) <- gsub(".canfam4", "", rownames(plink_in))

    rm(plink_data)

    res_ <- snp_bims[,c("chr", "pos", "alt", "ref")]
    colnames(res_) <- c("chr", "pos", "ALT", "REF")

    if(chr == "X"){
      res_$chr <- paste0("chrX")
      res_$chr_snp <- paste0("chrX_", res_$pos)
    } else {
      res_$chr <- paste0("chr",res_$chr)
      res_$chr_snp <- paste0(res_$chr, "_", res_$pos)
    }
    write.table(res_, paste0("metaData_SNPs/chr",chr,"_SNPs.tsv"))
    rm(res_)
  }
}

for (chr in c(1:38,"X")){
  print(paste0("chr",chr))
  for (j in batch){
    if(file.exists(paste0("SNP_metaData/chr",chr,"_SNP_batch",j,".tsv"))){next}
    SNPs <- read.table(paste0("metaData_SNPs/chr",chr,"_SNPs.tsv"))
    SNPs$chr_snp <- gsub("chrchr", "chr", SNPs$chr_snp)

    groups <- split(1:nrow(SNPs),
                    cut(1:nrow(SNPs), breaks = 1000, labels = FALSE)
    )

    # if(file.exists(paste0("SNP_metaData/SNP_batch",j,".tsv"))){next}

    SNPs <- SNPs[groups[[j]],]
    # SNPs <- SNPs[SNPs$region %in% c("chr3_42182260_42182338", "chr3_42237342_42237342", "chr3_42362631_42362768", "chr3_42385646_42385881" ,
    #                                                           "chr3_42386355_42386458", "chr3_42416665_42416668", "chr3_42465880_42466230", "chr3_42469976_42470076"),]

    regions <- na.omit(unique(SNPs[,c("chr_snp", "pos", "chr", "REF", "ALT"), drop = F]))
    colnames(regions) <- c("chr_snp", "snp", "chr", "REF", "ALT")
    regions$start <- as.numeric(regions$snp)
    regions$end<- as.numeric(regions$snp)

    chrs=unique(regions$chr)
    regions_list <- parallel::mclapply(chrs,
                                       function(x, regions){
                                         regions <- regions[regions$chr == x,]
                                         rownames(regions) <- paste0(regions$chr, "_", regions$snp, "_", regions$snp)
                                         return(regions)
                                       },
                                       regions = regions,
                                       mc.cores=n_cores)
    names(regions_list) <- chrs

    gtf_oi <- gtf[seqnames(gtf) %in% chrs,]
    transposon_annotation_oi <- transposon_annotation[transposon_annotation$seqnames %in% chrs,]

    SNP_metaData_TE_list <- parallel::mclapply(regions_list,
                                               DAPDNAm::TE_region_metaData_generation,
                                               transposons_annotation = transposon_annotation_oi,
                                               genome_gene_annotation = gtf_oi,
                                               sw_score_cutoff = 225,
                                               mc.cores=n_cores)

    for (i in 1:length(SNP_metaData_TE_list)) {
      if (i == 1){SNP_metaData_TE <- SNP_metaData_TE_list[[i]]} else {SNP_metaData_TE <- rbind(SNP_metaData_TE, SNP_metaData_TE_list[[i]])}
    }

    gtf_cpgisl_oi <- gtf_cpgisl[gtf_cpgisl$V2 %in% chrs,]

    SNP_metaData_list <- parallel::mclapply(regions_list,
                                            DAPDNAm::region_metaData_generation,
                                            genome_gene_annotation = gtf_oi,
                                            cpgisland_annotation = gtf_cpgisl_oi,
                                            mc.cores=n_cores)

    for (i in 1:length(SNP_metaData_list)) {
      if (i == 1){SNP_metaData <- SNP_metaData_list[[i]]} else {SNP_metaData <- rbind(SNP_metaData, SNP_metaData_list[[i]])}
    }

    chromatin_states_annotation_oi <- chromatin_states_annotation[seqnames(chromatin_states_annotation)%in% chrs,]

    SNP_metaData_ChromStates_list <- parallel::mclapply(regions_list,
                                                        DAPDNAm::chromatin_state_region_metaData_generation,
                                                        chromatin_states_annotation = chromatin_states_annotation_oi,
                                                        mc.cores=n_cores)

    for (i in 1:length(SNP_metaData_ChromStates_list)) {
      if (i == 1){SNP_metaData_ChromStates <- SNP_metaData_ChromStates_list[[i]]} else {SNP_metaData_ChromStates <- rbind(SNP_metaData_ChromStates, SNP_metaData_ChromStates_list[[i]])}
    }

    if (all(rownames(SNP_metaData_TE) == rownames(SNP_metaData))){SNP_metaData <- cbind(SNP_metaData, SNP_metaData_TE)} else {print("not all the transposable element rownames are the same!!!")}
    if (all(rownames(SNP_metaData_ChromStates) == rownames(SNP_metaData))){SNP_metaData <- cbind(SNP_metaData, SNP_metaData_ChromStates)} else {print("not all the chromatin states rownames are the same!!!")}
    # if (all(rownames(human_conserved_SNP_metaData) == rownames(SNP_metaData))){SNP_metaData <- cbind(SNP_metaData, human_conserved_SNP_metaData)} else {print("not all the human conserved regions rownames are the same!!!")}

    # Identify columns with the name "region"
    region_cols <- which(names(SNP_metaData) == "region")

    # Check if there are multiple "region" columns
    while (length(region_cols) > 1) {
      # Remove the second "region" column
      SNP_metaData <- SNP_metaData[, -region_cols[2]]
      # Identify columns with the name "region"
      region_cols <- which(names(SNP_metaData) == "region")
    }

    SNP_metaData$chr_snp <- gsub("^([^_]*_[^_]*)_.*$", "\\1", SNP_metaData$region)
    binding <- SNPs[SNPs$chr_snp %in% SNP_metaData$chr_snp,]
    binding <- binding[order(binding$chr_snp),]
    SNP_metaData <- SNP_metaData[order(SNP_metaData$chr_snp),]
    if(!all(binding$chr_snp == SNP_metaData$chr_snp)){stop("fix yo shit")}
    SNP_metaData<-cbind(SNP_metaData, binding[,c('chr','pos','ALT','REF')])

    # rownames(SNP_metaData) <- SNP_metaData$region

    all_results <- list()
    # load("inprogress.RData")

    library(VariantAnnotation)
    library(parallel)
    genome <- readDNAStringSet("/scratch/bmarine2/canFam4.fa", format = "fasta")

    chrs <- unique(SNP_metaData$chr)
    for (chr in chrs) {
      print(chr)

      chr_snps_oi <- unique(SNP_metaData$chr_snp[SNP_metaData$chr == chr])

      # Parallel processing using mclapply
      results <- mclapply(chr_snps_oi,
                          function(chr_snp, SNP_metaData, chr, genome) {
                            seq <- genome[[chr]]

                            oi <- SNP_metaData[SNP_metaData$chr_snp==chr_snp,,drop=F]
                            pos <- oi$pos
                            pos1 <- oi$pos+1
                            ref=oi$REF
                            alt=oi$ALT

                            posplus1REF <- as.character(subseq(seq, start = pos1, end = pos1))

                            pos1 <- oi$pos-1
                            posminus1REF <- as.character(subseq(seq, start = pos1, end = pos1))

                            motif <- as.character(subseq(seq, start = pos-50, end = pos+50))
                            motif_region <- paste0(chr, "_", pos-50, "_", pos+50)

                            return_ <- unique(data.frame(chr_snp = chr_snp, REF = ref[1], ALT = alt[1],
                                                         posplus1REF = posplus1REF[1],
                                                         posminus1REF = posminus1REF[1],
                                                         surr_motif = motif,
                                                         stringsAsFactors = FALSE))

                            return(return_)

                          },  SNP_metaData, chr, genome,
                          mc.cores = n_cores)

      # Combine results into a single data frame
      results <- do.call(rbind, results)

      # Deduplicate rows based on `chr_snp`
      results <- results[!duplicated(results$chr_snp), ]

      # Append to the list of results
      all_results[[chr]] <- results
    }

    if(length(chrs) == 1) {combined_results <- all_results[[1]]} else {
      # Combine all results from all chromosomes
      combined_results <- do.call(rbind, all_results)
    }

    # Remove duplicate rows if necessary
    combined_results <- combined_results[!duplicated(combined_results$chr_snp), ]

    # Merge combined results back into SNP_metaData
    SNP_metaData <- merge(SNP_metaData, combined_results[,c("chr_snp", "posplus1REF", "posminus1REF", "surr_motif")], by = "chr_snp", all.x = TRUE)

    CT_SNPs <- read.table("/scratch/bmarine2/CpG_gain_aware_fasta-250103/edited_snps.txt", header = F)

    SNP_metaData$CpG_disrupting <- "No"
    SNP_metaData$CpG_disrupting[SNP_metaData$chr_region %in% CT_SNPs$V1] <- "Yes"

    SNP_metaData$C_to_T <- "No"
    SNP_metaData$C_to_T[(SNP_metaData$REF == "C" & SNP_metaData$ALT == "T" | SNP_metaData$REF == "G" & SNP_metaData$ALT == "A")
    ] <- "Yes"

    write.table(x = SNP_metaData, file = paste0("SNP_metaData/chr",chr,"_SNP_batch",j,".tsv"))
  }
}



