#### Blaise Mariner now has a slight grasp on Sol and the dap mDNA data so now it is time to make a full workflow pipeline
#### this is part 2 of this series where we can visualize the regions, their coverage and methylation that we made from pt 1. 
#### part 3 will filter out the samples with a low lib size
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-01-12
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project

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
                      axis.text = element_text(family = "sans", size = 14, color="black"),axis.title.y = element_markdown(family = "sans", size = 14),   
                      axis.title.x = element_markdown(family = "sans", size = 14),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

#1.5 hrs ish with 25 cores (1 managing) and 250 Gb mem (200 errors out) during peak hours
cores_ <- 24

# this sets the working directory to this script's path
# if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# }
print(getwd())

# metadata
metaData <- read_rds("../metadata_samples/P1-DAP-metaData-2400216.rds")
p1 <- metaData
p1$weight <- as.numeric(p1$weight)
hist_plots <- list()
## age histograms by sex
p1$Age_at_sample <- as.numeric(p1$Age_at_sample)
# hist_plots[['age_female']] <-  ggplot(p1[p1$Sex == "Female",], aes(x = Age_at_sample)) +
#   geom_histogram(fill = "blue", color = "white") +
#   xlab("Age") +
#   labs(tag='A') +
#   ylab("Count") +
#   theme_blaise+
#   # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
#   scale_y_continuous(expand = c(0.01,0))+
#   scale_x_continuous(expand = c(0.005,0), limits = c(0,max(p1$Age_at_sample)+1))
# hist_plots[['age_male']] <- ggplot(p1[p1$Sex == "Male",], aes(x = Age_at_sample)) +
#   geom_histogram(fill = "red", color = "white") +
#   xlab("Age") +
#   ylab("Count") +
#   labs(tag='B') +
#   theme_blaise+
#   # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
#   scale_y_continuous(expand = c(0.01,0))+
#   scale_x_continuous(expand = c(0.005,0), limits = c(0,max(p1$Age_at_sample)+1))
# hist_male_age

## age histograms by sex and intact
p1$Age_at_sample <- as.numeric(p1$Age_at_sample)

hist_plots[['age_female']] <-ggplot(p1[p1$Sex == "Female" ,], aes(x = Age_at_sample)) +
  geom_histogram(fill = "blue", color = "white", bins = 40) +
  xlab("Age") +
  ylab("Count") +
  labs(tag='A') +
  ggtitle("")+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0), limits = c(0,50))+
  scale_x_continuous(expand = c(0.005,0), limits = c(0,max(p1$Age_at_sample)+1))

hist_plots[['age_male']] <- ggplot(p1[p1$Sex == "Male",], aes(x = Age_at_sample)) +
  geom_histogram(fill = "red", color = "white", bins = 40) +
  xlab("Age") +
  ylab("Count") +
  labs(tag='B') +
  ggtitle("")+
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0), limits = c(0,50))+
  scale_x_continuous(expand = c(0.005,0), limits = c(0,max(p1$Age_at_sample)+1))

hist_plots[['weight_female']] <- ggplot(p1[p1$Sex == "Female",], aes(x = weight)) +
  geom_histogram(fill = "blue", color = "white", bins = 20) +
  xlab("Weight") +
  ylab("Count") +
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0), limits = c(0,80))+
  scale_x_continuous(expand = c(0.005,0), limits = c(0,max(p1$weight, na.rm = T)))
# hist_plots[['weight_female']]

hist_plots[['weight_male']] <- ggplot(p1[p1$Sex == "Male",], aes(x = weight)) +
  geom_histogram(fill = "red", color = "white", bins = 20) +
  xlab("Weight") +
  ylab("Count") +
  theme_blaise+
  # geom_vline(xintercept = 5000, color = "red", linetype = "dashed", linewidth = 1.25) +
  scale_y_continuous(expand = c(0.01,0), limits = c(0,80))+
  scale_x_continuous(expand = c(0.005,0), limits = c(0,max(p1$weight, na.rm = T)))
# hist_plots[['weight_male']]

hist_plots[['weight_age_female']]  <- ggplot(p1[p1$Sex == "Female",], aes(x = Age_at_sample, y = weight)) +
  geom_point(color = "blue",size = 2) +
  xlab("Age") +
  ylab("Weight") +
  theme_blaise+
  scale_y_continuous(expand = c(0,0), limits = c(0,max(p1$weight, na.rm = T)))+
  scale_x_continuous(expand = c(0,0), limits = c(0,max(p1$Age_at_sample)+1))+
  geom_smooth(method = "loess", formula = y ~ x, color = "black")+
  coord_cartesian(clip = "off")

hist_plots[['weight_age_male']]  <- ggplot(p1[p1$Sex == "Male",], aes(x = Age_at_sample, y = weight)) +
  geom_point(color = "red",size = 2) +
  xlab("Age") +
  ylab("Weight") +
  theme_blaise+
  scale_y_continuous(expand = c(0,0), limits = c(0,max(p1$weight, na.rm = T)))+
  scale_x_continuous(expand = c(0,0), limits = c(0,max(p1$Age_at_sample)+1))+
  geom_smooth(method = "loess", formula = y ~ x, color = "black") +
  coord_cartesian(clip = "off")


plot1 <- (patchwork::wrap_plots(hist_plots, 
                            nrow = 3, 
                            ncol = 2, byrow = T))

# plot1
# head(metaData)
# devtools::install_github("blaisemariner17/DAPDNAm", force = TRUE)
library(DAPDNAm)

maxgap = 250
coverage_all_chr <- read_rds(paste0("../getcoverage_data_and_plots/maxGap_", maxgap, "_all_coverage_regions_oi.rds"))
methylation_all_chr <- read_rds(paste0("../getcoverage_data_and_plots/maxGap_", maxgap, "_all_methylation_regions_oi.rds"))

region_metaData <- readRDS("../metadata_regions/metaData_regions_maxgap_250.rds")
col_oi = c(
  "Promoter","exon", "intron",# "gene_bool",
  "CpG_island", "CpG_shore", "CpG_shelf",
  # "ChrSt_promoter", 
  "ChrSt_quies","ChrSt_heterochromatin",
  # "ChrSt_polycomb", 
  "ChrSt_enhancer",
  # 
  "TE"
)
percent_methylation_plots <- region_metaData_median_percent_methylation(coverage_all_chr, methylation_all_chr, region_metaData, 
                                                                        col_oi = col_oi, theme_blaise = theme_blaise
)

percent_methylation_plots[[1]] <-percent_methylation_plots[[1]] + labs(tag='C') 
plot2 <- (patchwork::wrap_plots(percent_methylation_plots, 
                            nrow = length(percent_methylation_plots), 
                            ncol = 1))

# plot1 | plot2

svglite("DAP_FIGURE1.svg", fix_text_size = F, width = 16)
print(plot1 | plot2)
dev.off()

