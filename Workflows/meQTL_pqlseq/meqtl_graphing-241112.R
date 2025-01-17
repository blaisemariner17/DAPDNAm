#!/usr/bin/env /packages/apps/spack/18/opt/spack/gcc-11.2.0/r-4.2.2-kpl/bin/Rscript
#### Blaise Mariner 
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-01-12
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project

rm(list = ls())  # Clear the workspace

# Determine if running in RStudio; set the number of cores accordingly
if (isRStudio <- Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Set working directory
  n_cores = 4  # Single-core for RStudio
  # load("in_progress.RData")
} else {
  n_cores = 20  # Multi-core for non-RStudio environments
}

# List of required libraries
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
  "patchwork",
  "ggrepel"
)
lapply(library_list, require, character.only = TRUE)  # Load required libraries

# Custom ggplot2 theme for consistency in plotting
theme_blaise <- theme(
  plot.title.position = "plot", 
  axis.text.x = element_text(angle = 0), 
  plot.title = element_markdown(family = "sans", size = 13, hjust = 0.5, color = "black", face = 'bold'), 
  plot.subtitle = element_text(family = "sans", size = 11, color = "black"), 
  axis.text = element_markdown(family = "sans", size = 14, color = "black"),
  axis.title.y = element_markdown(family = "sans", size = 20), 
  axis.title.x = element_markdown(family = "sans", size = 20), 
  panel.border = element_blank(), 
  axis.line = element_line(colour = "black", linewidth = 1), 
  axis.ticks = element_line(colour = "black", linewidth = 1), 
  legend.key.size = unit(1.5, 'cm'), 
  legend.key = element_rect(fill = NA), 
  legend.text = element_text(family = "sans", size = 20), 
  legend.title = element_blank(), 
  legend.background = element_blank(), 
  legend.box.background = element_blank(), 
  legend.text.align = 0, 
  panel.background = element_blank(), 
  panel.grid.major = element_line(colour = "black"), 
  panel.grid.minor = element_blank()
) + removeGrid()

# Load metadata for samples and filter based on use criteria
metaData <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/metaData_samples/all_precision_metaData.rds")
metaData <- metaData[metaData$use_for_baseline_precision_clock == "yes", ]
metaData <- metaData[order(metaData$dog_id), ]

# Load PQLseq results and format the dataframe
pqlseq_res <- read.table("all_pqlseq_res.tsv")
print(dim(pqlseq_res))
pqlseq_res <- unique(pqlseq_res)
colnames(pqlseq_res) <- c(
  "idfk", "region_snp", "nSNP" ,"interceptSNP" ,"se_interceptSNP", "betaSNP" ,
  "se_betaSNP", "pvalueSNP" ,"h2SNP" ,"sigma2SNP", "convergedSNP", "elapsed_timeSNP"
)
rownames(pqlseq_res) <- pqlseq_res$region_snp

dim(pqlseq_res)
dim(pqlseq_res[pqlseq_res$convergedSNP == T,])
head(pqlseq_res)

# Filter results to include only converged models
for (colname in colnames(pqlseq_res)[grepl("converged", colnames(pqlseq_res))]) {
  print(colname)
  print("FALSE")
  print((dim(pqlseq_res[pqlseq_res[, c(paste(colname))] == F, ])))
  print("TRUE")
  print((dim(pqlseq_res[pqlseq_res[, c(paste(colname))] == T, ])))
  pqlseq_res <- pqlseq_res[pqlseq_res[, c(paste(colname))] == TRUE, ]
}

# Parse regions and SNPs from `region_snp` column
tmp <- as.data.frame(do.call("rbind", stringr::str_split(pqlseq_res$region_snp, "_")))
colnames(tmp) <- c("chr", "start", "end", "snp")
tmp$region <- paste(tmp$chr, tmp$start, tmp$end, sep = "_")
tmp$start <- as.numeric(tmp$start)
tmp$end <- as.numeric(tmp$end)
tmp$region_length <- tmp$end - tmp$start
tmp$chr_snp <- paste0(tmp$chr, "_", tmp$snp)
tmp$regions_snps <- paste(tmp$region, tmp$snp, sep = "_")
pqlseq_res <- cbind(pqlseq_res, tmp)
rm(tmp)
head(pqlseq_res)

# Unique regions and chr_snp identifiers
unique_regions <- unique(pqlseq_res$region)
length(unique_regions)
unique_chrsnps<- unique(pqlseq_res$chr_snp)
length(unique_chrsnps)

library(parallel)  # Load parallel processing library

# Function to compute SNPs per region
compute_nsnps <- function(region, regions_snps) {
  nsnps <- nrow(regions_snps[regions_snps$region == region, ])
  length <- unique(regions_snps$region_length[regions_snps$region == region])
  return(data.frame("region" = region, "nsnps" = nsnps, "region_length" = length))
}

# Parallel processing to calculate SNPs per region
nsnps_per_region <- mclapply(unique_regions, compute_nsnps, regions_snps = pqlseq_res, mc.cores = 20)
nsnps_per_region <- do.call(rbind, nsnps_per_region)

# Function to calculate mean effects per SNP
compute_snps_info <- function(snp, pqlseq_res) {
  mean_effect_snp <- mean(pqlseq_res$betaSNP[pqlseq_res$chr_snp == snp])
  mean_effect_age <- mean(pqlseq_res$betaAge_at_sample[pqlseq_res$chr_snp == snp])
  return(data.frame("snp" = snp, "mean_effect_snp" = mean_effect_snp, "mean_effect_age" = mean_effect_age))
}

# Parallel processing to calculate mean effects for SNPs
snps_info <- mclapply(unique_chrsnps, compute_snps_info, pqlseq_res = pqlseq_res, mc.cores = 20)
snps_info <- do.call(rbind, snps_info)

# Data for ggplot visualization
for_ggplot <- nsnps_per_region

# Plot 1: Histogram of SNPs associated per region
plot1 <- ggplot(for_ggplot, aes(x = nsnps)) +
  geom_histogram(color = "white", fill = "#736372") +
  theme_blaise +
  xlab("No. SNPs associated per region") +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0))

# Plot 2: Histogram of regions within 50kb of each SNP
for_ggplot <- data.frame(table(pqlseq_res$chr_snp))
plot2 <- ggplot(for_ggplot, aes(x = Freq)) +
  geom_histogram(color = "white", fill = "#AAC0AA") +
  theme_blaise +
  xlab("No. regions within 50kb of each SNP") +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0))

# Save combined plots 1 and 2 as an SVG
svglite("meQTLplotting/plots12.svg", fix_text_size = FALSE)
print(plot1 + plot2)
dev.off()

plot3 <- ggplot(pqlseq_res, aes(x = betaSNP)) +
  theme_blaise + 
  geom_histogram(color = "white", fill = "#7A918D") + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkred", linewidth = 2)+
  xlab("β<sub>SNP</sub>") + ylab("Count") +
  scale_y_continuous(expand = c(0,0))

plot4 <- ggplot(pqlseq_res, aes(x = pvalueSNP)) +
  theme_blaise + 
  geom_histogram(color = "white", fill = "#7A918D") + 
  xlab("pval<sub>SNP</sub>") + ylab("Count") +
  scale_y_continuous(expand = c(0,0))

svglite("meQTLplotting/plots34.svg", fix_text_size = F)
print(plot3 + plot4)
dev.off()

######### What about the distances between SNP site and region-- what can that tell us about the relationship
pqlseq_res$snp<- as.numeric(pqlseq_res$snp)
pqlseq_res$reg_start_snp <- abs(pqlseq_res$start - pqlseq_res$snp)
pqlseq_res$reg_end_snp <- abs(pqlseq_res$end - pqlseq_res$snp)
pqlseq_res$distance_reg_snp <- pqlseq_res$reg_start_snp
pqlseq_res$distance_reg_snp[pqlseq_res$reg_end_snp <= pqlseq_res$reg_start_snp] <- pqlseq_res$reg_end_snp[pqlseq_res$reg_end_snp <= pqlseq_res$reg_start_snp]
pqlseq_res$distance_reg_snp[pqlseq_res$snp > pqlseq_res$start & pqlseq_res$snp < pqlseq_res$end] <- 0

# number of SNPs within a region: 
sum(pqlseq_res$distance_reg_snp == 0)

#this will be done from running the sbatch script next time.

pqlseq_res$padjSNP <- p.adjust(pqlseq_res$pvalueSNP, method = "BH")

snp_sig_regions <- unique(pqlseq_res$region[pqlseq_res$padjSNP < 0.05])
pqlseq_res$sig_snp <- "no"
pqlseq_res$sig_snp[pqlseq_res$region %in% snp_sig_regions] <- "yes"

for_ggplot <- pqlseq_res

plot6 <- ggplot(for_ggplot[
  sample(x = 1:nrow(for_ggplot), 10000),], 
  aes(x = distance_reg_snp)) +
  theme_blaise + 
  geom_density()+ 
  xlab("distance b/w SNP and region (bp)") + ylab("Density") +
  scale_y_continuous(expand = c(0,0))

plot7 <- ggplot(for_ggplot, 
                aes(y = -log10(pvalueSNP), x = distance_reg_snp)) +
  theme_blaise + 
  geom_point(alpha = 0.5)+ 
  xlab("distance b/w SNP and region (bp)") + ylab("-log<sub>10</sub>(p-value)") +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off")

svglite("meQTLplotting/plots67.svg", fix_text_size = F)
print(plot6 + plot7)
dev.off()

print(paste0(nrow(for_ggplot[for_ggplot$distance_reg_snp == 0,]), " SNPs are within a region."))

summary(lm(-log10(pvalueSNP) ~ distance_reg_snp, data = for_ggplot))

testing_distance <- for_ggplot#[for_ggplot$sig_snp == 'yes',]

plot_test <- ggplot(testing_distance,
       aes(x = pvalueSNP)) +
  theme_blaise +
  # ggtitle("significant SNP-region pairs only") +
  geom_histogram(color = "white", fill = "black")+
  xlab("pval<sub>SNP</sub>") + ylab("Count") 

plot_test_1 <- ggplot(testing_distance[testing_distance$distance_reg_snp < 10000,], 
                      aes(x = pvalueSNP)) +
  theme_blaise + 
  ggtitle("SNP-region pairs <10kb apart") +
  geom_histogram(color = "white", fill = "black")+ 
  xlab("pval<sub>SNP</sub>") + ylab("Count") 
  
plot_test_2 <- ggplot(testing_distance[testing_distance$distance_reg_snp > 10000 & testing_distance$distance_reg_snp < 25000,], 
                      aes(x = pvalueSNP)) +
  theme_blaise + 
  ggtitle("SNP-region pairs 10kb< <25kb apart") +
  geom_histogram(color = "white", fill = "black")+ 
  xlab("pval<sub>SNP</sub>") + ylab("Count") 

plot_test_3 <- ggplot(testing_distance[testing_distance$distance_reg_snp > 25000,], 
                      aes(x = pvalueSNP)) +
  theme_blaise + 
  ggtitle("SNP-region pairs >25kb apart") +
  geom_histogram(color = "white", fill = "black")+ 
  xlab("pval<sub>SNP</sub>") + ylab("Count") 


plot_test

plot_test_1 | plot_test_2 | plot_test_3


#######################################################
#######################################################
#######################################################

og_bmm <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/pqlseq_results/pqlseq_res_1n_maxgap250.rds")
for (colname in colnames(og_bmm)[grepl("converged", colnames(og_bmm))]){
  og_bmm <- og_bmm[og_bmm[,c(paste(colname))] == T,]
}
og_bmm_sig <- rownames(og_bmm)[og_bmm$padj_Age_at_sample < 0.05]
age_sig_regions <- rownames(og_bmm)[og_bmm$padj_Age_at_sample < 0.05]

all_snps <- unique(pqlseq_res$chr_snp)
sig_snps <- unique(pqlseq_res$chr_snp[pqlseq_res$padjSNP < 0.05])
notsig_snps <- unique(pqlseq_res$chr_snp[! pqlseq_res$chr_snp %in% sig_snps])
sig_snps_age <- unique(pqlseq_res$chr_snp[pqlseq_res$region %in% age_sig_regions])


a <- sig_snps[sig_snps %in% sig_snps_age]
b <- sig_snps[!sig_snps %in% sig_snps_age]
c <- notsig_snps[notsig_snps %in% sig_snps_age]
d <- all_snps[!all_snps %in% c(a,b,c)]
 
data_oi <- matrix(c(
  length(a),
  length(b),
  length(c),
  length(d)
  ),
                  nrow=2, ncol=2, byrow=TRUE)
colnames(data_oi) <- c("SNP_that_influences_meth", "not")
rownames(data_oi) <- c("Region_within_50kb_of_SNP_ass_with_age", "not")
data_oi

odds_ratio <- epitools::oddsratio(data_oi)
data.frame("class" = "SNPs_assc_with_meth_AND_age_related_meth", 
           "pval" = odds_ratio$p.value[2,2],
           "odds_ratio" = odds_ratio$measure[2,1],
           "lower" = (odds_ratio$measure[2,2]), 
           "upper" = (odds_ratio$measure[2,3]))

###### Let's look a little more at the results
save.image("in_progress.RData")

col_oi = c(
  "Promoter", #"gene_bool",
  "exon", "intron",
  "CpG_island", "CpG_shore", "CpG_shelf",
  # "ChrSt_promoter",
  "ChrSt_quies","ChrSt_heterochromatin",
  "ChrSt_polycomb",
  "ChrSt_enhancer",
  "DNA_transposon", "Retrotransposon",
  "TE", "LINE1"
)

######## we want to know where these SNPs are-- like where are they located and everything...

snp_metaData <- readRDS(paste0("metaData_SNPs/metaData_SNPs.rds"))
snp_metaData$snp <- gsub(".*_", "", snp_metaData$region)
snp_metaData$chr <- gsub("_.*", "", snp_metaData$region)
snp_metaData$chr_snp <- paste0(snp_metaData$chr, "_", snp_metaData$snp)
snp_metaData$fdr_snp <- "not significant"
snp_metaData$fdr_snp[snp_metaData$chr_snp %in% unique(pqlseq_res$chr_snp[pqlseq_res$sig_snp == "yes"])] <- "SNP influence"

table(snp_metaData$fdr_snp)

snp_metaData <- snp_metaData[!rowSums(snp_metaData[,paste(col_oi)])==0,]
i=1
for (col in c(col_oi)){
  a <- nrow(snp_metaData[snp_metaData[,paste(col)] == 1  & snp_metaData$fdr_snp == "SNP influence",])
  b <- nrow(snp_metaData[snp_metaData[,paste(col)] == 1 & snp_metaData$fdr_snp == "not significant",])
  c <- nrow(snp_metaData[(!snp_metaData[,paste(col)] == 1)  & snp_metaData$fdr_snp == "SNP influence",])
  d <- nrow(snp_metaData[(!snp_metaData[,paste(col)] == 1) & snp_metaData$fdr_snp == "not significant",])
  
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

############################################################################################################
cov <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_coverage_regions_1n_coverage.rds")
meth <- readRDS("/scratch/bmarine2/BASELINE-PRECISION-240820/getcoverage_data_and_plots/maxGap_250_all_methylation_regions_1n_coverage.rds")
perc_meth <- meth / cov
perc_meth <- perc_meth[,colnames(perc_meth) %in% rownames(metaData)]
perc_meth <- perc_meth[,order(colnames(perc_meth))]
metaData <- metaData[rownames(metaData) %in% colnames(perc_meth),]
metaData <- metaData[order(rownames(metaData)),]
pqlseq_res$sig_age<- "no"
pqlseq_res$sig_age[pqlseq_res$region %in% age_sig_regions] <- "yes"
oi <- pqlseq_res[pqlseq_res$chr_snp %in% sig_snps,]

oi <- oi[order(oi$pvalueSNP),]
head(oi)

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

for (i in c(10:10)){
  chr <- oi$chr[i]
  chr_snp <- oi$chr_snp[i]
  region <- oi$region[i]
  rowname <- rownames(oi)[i]
  chr_region_snp <- oi$region_snp[i]
  age_ass <- oi$sig_age[i]
  
  if(age_ass == "yes"){age_ass <- "*"}
  if(age_ass == "no"){age_ass <- ""}
  
  perc_meth_oi <- perc_meth[region,, drop = F]

  metaData_oi <- cbind(metaData, t(perc_meth_oi))
  colnames(metaData_oi) <- c(colnames(metaData), "region_oi")
  
  plink_in <- MultiPhen::read.plink(paste0("/scratch/bmarine2/meQTL-241101-GEMMA/gemma_out/filtered_genotypes_chr", gsub("chr","",chr)))
  plink_in <- (plink_in[,colSums(plink_in == 2) < .9*nrow(plink_in) & colSums(plink_in == 1) < .9*nrow(plink_in)])
  colnames(plink_in) <- paste0("chr", "", colnames(plink_in))
  plink_in <- plink_in[order(rownames(plink_in)),]
  plink_in <- plink_in[rownames(plink_in) %in% metaData_oi$dog_id,]
  
  metaData_oi <- metaData_oi[metaData_oi$dog_id %in% rownames(plink_in),]
  metaData_oi <- metaData_oi[order(metaData_oi$dog_id),]
  metaData_oi <- metaData_oi[metaData_oi$dog_id %in% rownames(plink_in),]
  metaData_oi <- metaData_oi[order(metaData_oi$dog_id),]
  metaData_oi$SNP <- (plink_in[,c(paste(chr_snp))])
  metaData_oi <- metaData_oi[metaData_oi$SNP %in% c(0,1,2),]
  
  max_age <- max(metaData_oi$Age_at_sample)
  
  plot1 <- ggplot(metaData_oi[metaData_oi$Sex == "Male",], 
                  aes(color = SNP, y = region_oi, x = factor(SNP, levels = c(0,1,2)))) +
    ggbeeswarm::geom_beeswarm() +
    geom_boxplot(width = 0.2, outliers = F, color = "black")+
    theme_blaise + theme(legend.position = "none") +
    ggtitle(paste0("Male; ", chr_snp)) +
    ylab(paste0(region, age_ass)) +
    xlab(paste0(chr_snp)) +
    geom_smooth(method = "lm")+
    scale_y_continuous(limits = c(0,1), labels = scales::percent, expand = c(0,0))+
    # scale_x_continuous(limits = c(0, max_age), expand = c(0,0))+
    coord_cartesian(clip = "off")
  
  plot2 <- ggplot(metaData_oi[metaData_oi$Sex == "Female",], 
                  aes(color = SNP, y = region_oi, x = factor(SNP, levels = c(0,1,2)))) +
    ggbeeswarm::geom_beeswarm() +
    geom_boxplot(width = 0.2, outliers = F, color = "black")+
    theme_blaise + theme(legend.position = "none") +
    ggtitle(paste0("Female; ", chr_snp)) +
    ylab(paste0(region, age_ass)) +
    xlab(paste0(chr_snp)) +
    geom_smooth(method = "lm")+
    scale_y_continuous(limits = c(0,1), labels = scales::percent, expand = c(0,0))+
    # scale_x_continuous(limits = c(0, max_age), expand = c(0,0)) +
    coord_cartesian(clip = "off")
  
  print(plot2 + plot1)
}




