#### Blaise Mariner now has a slight grasp on Sol and the dap mDNA data so now it is time to make a full workflow pipeline
#### this is part 2 of this series where we can visualize the regions, their coverage and methylation that we made from pt 1.
#### part 3 will filter out the samples with a low lib size
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2023-11-09
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project
rm(list=ls())

library_list <- c(
  "ggplot2",
  "svglite",
  "ggExtra",
  "ggtext",
  "tidyverse",
  "corrplot"
)
lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 24, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      axis.text = element_text(family = "sans", size = 18, color="black"),       axis.title = element_text(family = "sans", size = 20, color="black"),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

cores_ <- 1

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

# metadata
p1 <- read_rds("metadata_samples/P1-DAP-metaData-240510.rds")
p2 <- read_rds("metadata_samples/P2-DAP-metaData-240510.rds")
p3 <- read_rds("metadata_samples/P3-DAP-metaData-240510.rds")

metaData <- rbind(p1,p2,p3);nrow(metaData)
# head(metaData)

#get the maxgaps I used from the previous loop
maxgaps_split <- str_split(list.files(path = "getcoverage_data_and_plots/", pattern = "_all_coverage_regions_oi.rds"), fixed("_"))
maxgaps <- c()
for (i in 1:length(maxgaps_split)) {
  maxgap <- maxgaps_split[[i]][2]
  maxgaps <- append(maxgaps, maxgap)
}
maxgaps = 250

# define the number of cores
cores_ <- min(length(maxgaps),cores_)



#read in the data
regions_all_chr <- readRDS(paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_regions_oi.rds"))
coverage_all_chr <- readRDS(paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_coverage_regions_oi.rds"))
methylation_all_chr <- readRDS(paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_methylation_regions_oi.rds"))

if (dir.exists(paste0("getcoverage_data_and_plots/", maxgap)) == FALSE) {
  dir.create(paste0("getcoverage_data_and_plots/", maxgap))
}

# region length histogram
# plot the length of each each region as a histogram
regions_hist_length <- ggplot(regions_all_chr[regions_all_chr$length < 3000,], aes(x = length)) +
  geom_histogram(fill = "black", color = "white") +
  xlab("length") +
  ylab("Count") +
  ggtitle(paste0("Region length maxgap = ", maxgap))+
  labs(caption = (paste0("There are ", nrow(regions_all_chr[regions_all_chr$length > 3000,]), " regions longer than 3,000 bp")))+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap,"/hist_region_length.svg"), fix_text_size = FALSE)
print(regions_hist_length)
dev.off()

# plot the n cpgs in each region as a histogram
regions_hist_nCpG <- ggplot(regions_all_chr[regions_all_chr$n<100,], aes(x = n)) +
  geom_histogram(fill = "black", color = "white") +
  xlab("Number of CpGs by region") +
  ylab("Count") +
  ggtitle(paste0("Number CpGs by region maxgap = ", maxgap))+
  labs(caption = (paste0("There are ", nrow(regions_all_chr), " total regions and ", nrow(regions_all_chr[regions_all_chr$n>100,])," regions with more than 100 nCpGs")))+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/hist_nCpG_by_region.svg"), fix_text_size = FALSE)
print(regions_hist_nCpG)
dev.off()

regions_hist_length_nCpG <- ggplot(regions_all_chr, aes(x = length, y = n)) +
  geom_point(fill = "black", color = "black") +
  xlab("region length") +
  ylab("number of CpGs in the region") +
  ggtitle(paste0("maxgap = ", maxgap))+
  # labs(caption = (paste0("There are ", nrow(regions_all_chr[regions_all_chr$length > 100000,]), " regions longer than 100,000 bp")))+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/region_length_nCpG.svg"), fix_text_size = FALSE)
print(regions_hist_length_nCpG)
dev.off()

regions_all_chr$length_by_nCpG <- regions_all_chr$length / regions_all_chr$n

hist_length_by_nCpG <- ggplot(regions_all_chr, aes(x = length_by_nCpG)) +
  geom_histogram(fill = "black", color = "white") +
  xlab("Length / number of CpGs by region") +
  ylab("Count") +
  ggtitle(paste0("Length / number CpGs by region"))+
  # labs(caption = (paste0("There are ", nrow(regions_all_chr[regions_all_chr$n>500,])," regions with more than 500 nCpGs")))+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/hist_length_by_nCpG_by_region.svg"), fix_text_size = FALSE)
print(hist_length_by_nCpG)
dev.off()

#### we can also look at some perc methylation

perc_methylation <- as.data.frame(rowSums(methylation_all_chr) / rowSums(coverage_all_chr))
perc_methylation$pm_byregion <- perc_methylation[,1]

hist_perc_meth_by_region <- ggplot(perc_methylation, aes(x = pm_byregion)) +
  geom_histogram(fill = "black", color = "white") +
  xlab("Mean percent methylation by region") +
  ylab("Count") +
  ggtitle(paste0("Mean percent methylation by region"))+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/hist_perc_meth_by_region.svg"), fix_text_size = FALSE)
print(hist_perc_meth_by_region)
dev.off()

#### Let's look at some lid histograms
coverage <- as.data.frame(colSums(coverage_all_chr))
colnames(coverage) <- c("coverage")

hist_lid_lib_size <- ggplot(coverage, aes(x = coverage)) +
  geom_histogram(fill = "black", color = "white") +
  xlab("Total coverage by lid sample") +
  ylab("Count") +
  ggtitle(paste0("Total coverage by lid sample"))+
  labs(caption = (paste0("Line = 1.5x10<sup>7</sup> coverage")))+
  theme_blaise+
  # geom_vline(xintercept = 150, color = "red", linetype = "dashed", linewidth = 1.25) +
  # geom_vline(xintercept = 15000000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/hist_lid_sample_coverage.svg"), fix_text_size = FALSE)
print(hist_lid_lib_size)
dev.off()

# and perc_methylation by sample

# perc_methylation <- as.data.frame(colSums(methylation_all_chr) / colSums(coverage_all_chr))
# perc_methylation$pm_bysample <- perc_methylation[,1]
#
# hist_perc_meth_by_sample <- ggplot(perc_methylation, aes(x = pm_bysample)) +
#   geom_histogram(fill = "black", color = "white") +
#   xlab("Mean percent methylation by lid sample") +
#   ylab("Count") +
#   ggtitle(paste0("Mean percent methylation by lid sample; maxgap = ", maxgap))+
#   theme_blaise+
#   # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
#   scale_y_continuous(expand = c(0.01,0))+
#   scale_x_continuous(expand = c(0.005,0), limits = c(0,1))
#
# svglite(paste0("getcoverage_data_and_plots/", maxgap, "/hist_perc_meth_by_sample.svg"), fix_text_size = FALSE)
# print(hist_perc_meth_by_sample)
# dev.off()

perc_methylation <- as.data.frame(rowSums(methylation_all_chr) / rowSums(coverage_all_chr))
perc_methylation$pm_byregion <- perc_methylation[,1]
all(rownames(perc_methylation) == rownames(regions_all_chr))
perc_methylation$region_length <- regions_all_chr$length
perc_methylation <- as.data.frame(perc_methylation)

perc_meth_by_sample_w_regionsize <- ggplot(perc_methylation, aes(y = region_length, x = pm_byregion)) +
  geom_point(fill = "black", color = "black", alpha=0.1) +
  ylab("region length") +
  xlab("percent CpG's methylated") +
  ggtitle(paste0(""))+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_log10(expand = c(0.01,0))+
  scale_x_continuous(expand = c(0.005,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/mean_perc_meth_by_sample_div_regionsize.svg"), fix_text_size = FALSE)
print(perc_meth_by_sample_w_regionsize)
dev.off()

### generating a correlation matrix of percent methylation
# and then run a pca on the correlation matrix to determine the axes of variation
sampling <- sample(nrow(methylation_all_chr), 10000)

perc_methylation <- as.matrix(methylation_all_chr / coverage_all_chr)[sampling,]

perc_methylation <- perc_methylation[,colnames(perc_methylation) %in% metaData$lid_pid]
perc_methylation <- perc_methylation[,order(colnames(perc_methylation))]
metaData <- metaData[metaData$lid_pid %in% colnames(perc_methylation),]
metaData <- metaData[order(metaData$lid_pid),]
all(metaData$lid_pid == colnames(perc_methylation))

perc_meth_batch_rm <- limma::removeBatchEffect(perc_methylation,
                                               metaData$prep_date)

corr_mat <- stats::cor(perc_meth_batch_rm,use="pairwise.complete.obs")
write_rds(x = corr_mat, file = paste0("getcoverage_data_and_plots/", maxgap, "/correlation_matrix-20231206.rds"))

pca_corr_mat <- prcomp(corr_mat)

my_pca <- as.data.frame(pca_corr_mat$rotation)
metaData_ <- metaData[metaData$lid_pid %in% rownames(my_pca),]
my_pca <- my_pca[rownames(my_pca) %in% metaData_$lid_pid,]
metaData_ <- metaData_[order(metaData_$lid_pid),]

my_pca.var <- pca_corr_mat$sdev ^ 2
propve <- my_pca.var / sum(my_pca.var)

for_ggplot <- as.data.frame(propve[1:5])
for_ggplot$index <- 1:nrow(for_ggplot)
colnames(for_ggplot) <- c("prop_var", "PC")

scree_ <- ggplot(for_ggplot, aes(x = PC, y = prop_var)) +
  geom_col(fill = "black") +
  ggtitle(paste0("Scree plot; maxgap = ", maxgap))+
  ylab("Proportion of the variance explained")+
  xlab("Principal component")+
  theme_blaise+
  scale_y_continuous(expand = c(0.01,0), limits = c(0,1))+
  scale_x_continuous(expand = c(0.01,0))

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/scree_plot-corrmat.svg"), fix_text_size = F)
print(scree_)
dev.off()

my_pca$Age_at_sample <- as.numeric(metaData_$Age_at_sample)
my_pca$prep_date <- as.factor(metaData_$Prep_date)
my_pca$weight <- as.numeric(metaData_$weight)
my_pca$Breed <- as.factor(metaData_$Breed)
my_pca$Sex <- as.factor(metaData_$Sex)
my_pca$Breed_size_class <- as.factor(metaData_$Breed_Size_Class_at_HLES)
my_pca$Sex_class <- as.factor(metaData_$Sex_Class)
my_pca$Breed_status <- as.factor(metaData_$Breed_status)

plot_ <- ggplot(my_pca, aes(x = PC1, y = PC2)) +
  geom_point() +
  xlab(paste0("PC1: ", round(propve[1], 3) * 100, "% variance")) +
  ylab(paste0("PC2: ", round(propve[2], 3) * 100, "% variance")) +
  ggtitle(paste0(""))+
  theme_blaise + theme(legend.position = 'none')

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/PCA_corr_color-prep-date.svg"), fix_text_size = F)
print(plot_)
dev.off()

# plot_ <- ggplot(my_pca, aes(x = PC1, y = PC2, color = Age_at_sample)) +
#   geom_point() +
#   xlab(paste0("PC1: ", round(propve[1], 3) * 100, "% variance")) +
#   ylab(paste0("PC2: ", round(propve[2], 3) * 100, "% variance")) +
#   ggtitle(paste0("maxgap = ", maxgap))+
#   theme_blaise
#
# svglite(paste0("getcoverage_data_and_plots/", maxgap, "/PCA_corr_color-age-at-prep.svg"), fix_text_size = F)
# print(plot_)
# dev.off()
#
# plot_ <- ggplot(my_pca, aes(x = PC1, y = PC2, color = weight)) +
#   geom_point() +
#   xlab(paste0("PC1: ", round(propve[1], 3) * 100, "% variance")) +
#   ylab(paste0("PC2: ", round(propve[2], 3) * 100, "% variance")) +
#   ggtitle(paste0("maxgap = ", maxgap))+
#   theme_blaise
#
# svglite(paste0("getcoverage_data_and_plots/", maxgap, "/PCA_corr_color-weight.svg"), fix_text_size = F)
# print(plot_)
# dev.off()
#
# plot_ <- ggplot(my_pca, aes(x = PC1, y = PC2, color = Sex)) +
#   geom_point() +
#   xlab(paste0("PC1: ", round(propve[1], 3) * 100, "% variance")) +
#   ylab(paste0("PC2: ", round(propve[2], 3) * 100, "% variance")) +
#   ggtitle(paste0("maxgap = ", maxgap))+
#   theme_blaise
#
# svglite(paste0("getcoverage_data_and_plots/", maxgap, "/PCA_corr_color-sex.svg"), fix_text_size = F)
# print(plot_)
# dev.off()

metaData_$fixed_num[metaData_$fixed == "Spayed"] <- 1
metaData_$fixed_num[metaData_$fixed == "Neutered"] <- 1
metaData_$fixed_num[metaData_$fixed == "Intact"] <- 0
metaData_$fixed_num[metaData_$fixed == "unknown"] <- NA

metaData_$Sex[metaData_$Sex == "Male"] <- 0
metaData_$Sex[metaData_$Sex == "Female"] <- 1
metaData_$Sex[metaData_$Sex == "unknown"] <- NA

metaData_$Breedclass_num[metaData_$Breed_status == "unknown"] <- NA
metaData_$Breedclass_num[metaData_$Breed_status == "Mixed breed"] <- 0
metaData_$Breedclass_num[metaData_$Breed_status == "Purebred"] <- 1

for_corrplot <- as.data.frame(my_pca[,c("PC1", "PC2", "PC3")])
all(rownames(for_corrplot) == metaData_$lid_pid)
# for_corrplot$Unique_reads <- as.numeric(metaData_$Unique_reads)
for_corrplot$Perc_methylated <- as.numeric(metaData_$perc_meth)
for_corrplot$Perc_unique_reads <- as.numeric(metaData_$Perc_unique)
# for_corrplot$chrX_ratio <- as.numeric(metaData_$chrX_ratio)
for_corrplot$Age_at_sample <- as.numeric(metaData_$Age_at_sample)
for_corrplot$height_adjusted_size <- as.numeric(metaData_$height_adjusted_size)
for_corrplot$Sex <- metaData_$Sex_num
for_corrplot$Fixed <- metaData_$fixed_num
for_corrplot$Purebred <- metaData_$Breedclass_num
# for_corrplot$Prep_date <- metaData_$Prep_date

M <- cor(for_corrplot, use="pairwise.complete.obs")
testRes = cor.mtest(for_corrplot, conf.level = 0.95)

svglite(paste0("getcoverage_data_and_plots/", maxgap, "/correlationplot_hm.svg"), fix_text_size = F)
corrplot(M,
         # p.mat = testRes$p, sig.level = 0.05,
         type = "lower", method = 'color',
         order = "original",
         col = COL2('PRGn'),
         tl.col = 'black', tl.srt = 45, tl.pos = "l",
         # title = paste0("Correlation plot for maxgap = ", maxgap, ".")
         addrect = 3
)
dev.off()


