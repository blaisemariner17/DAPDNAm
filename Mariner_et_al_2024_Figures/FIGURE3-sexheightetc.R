#### Blaise Mariner now has a slight grasp on Sol and the dap mDNA data so now it is time to make a full workflow pipeline
#### this is part 2 of this series where we can visualize the regions, their coverage and methylation that we made from pt 1. 
#### part 3 will filter out the samples with a low lib size
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-01-12
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project
col_oi = c(
  "Promoter", "gene_bool",
  "CpG_island", "CpG_shore", "CpG_shelf",
  # "ChrSt_promoter", 
  "ChrSt_quies","ChrSt_heterochromatin",
  "ChrSt_polycomb", "ChrSt_enhancer",
  # "exon", "intron",
  "DNA.transposon", "Retrotransposon", "TE"
)
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
region_metaData <- readRDS("../metadata_regions/metaData_regions_maxgap_250.rds")

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
  "ggbeeswarm",
  "ggsignif"
)
lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 12, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      
                      axis.text = element_text(family = "sans", size = 14, color="black"),axis.title.y = element_markdown(family = "sans", size = 14),   
                      axis.title.x = element_markdown(family = "sans", size = 14),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

#1.5 hrs ish with 25 cores (1 managing) and 250 Gb mem (200 errors out) during peak hours
cores_ <- 24

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())
set.seed(100)

# metadata
metaData <- read_rds("../metadata_samples/P1-DAP-metaData-2400216.rds")

fit <- readRDS("../pqlseq_results/pqlseq_res_breedsizeintage_maxgap250.rds")
fit <- fit[!grepl("chrY", rownames(fit)),]
pqlseq_res <- readRDS("../pqlseq_results/pqlseq_res_maxgap250.rds")
pqlseq_res <- pqlseq_res[order(pqlseq_res$padj),]
pqlseq_res <- pqlseq_res[pqlseq_res$converged == T,]

pqlseq_res$count <- 1:nrow(pqlseq_res)
pqlseq_res$fdr <- pqlseq_res$count * pqlseq_res$padj
pqlseq_res$fdr_perc <- pqlseq_res$fdr / pqlseq_res$count
fit <- fit[rownames(fit) %in% rownames(pqlseq_res)[pqlseq_res$fdr_perc < 0.01 & pqlseq_res$beta<0],]
col_oi <- "TE"
for (col in c("Small", "Medium", "Standard","Large", "Giant")){
  
  fit_ <- fit[,grepl(col, colnames(fit))]
  
  new_colname <-  paste0("Age_at_sample_int_", gsub(".*size","",col))
  
  colnames(fit_) <- gsub(paste0("_",new_colname), "", colnames(fit_))
  
  fit_ <- fit_[order(fit_$padj),]
  fit_ <- fit_[fit_$converged == T,]
  
  fit_$count <- 1:nrow(fit_)
  fit_$fdr <- fit_$count * fit_$padj
  fit_$fdr_perc <- fit_$fdr / fit_$count
  
  age_effect_plots_fdr05=list()
  region_metaData_oi <- region_metaData[region_metaData[,c(paste(col_oi))] == 1,]
  fit_oi <- fit_[rownames(fit_) %in% region_metaData_oi$region,]
  
  fit_oi$group <- col
  
  if (col == "Small") { fit_oi_all <- fit_oi} else {fit_oi_all <- rbind(fit_oi_all, fit_oi)}
}

# fit_oi_all$beta <- abs(fit_oi_all$beta)
sampling <- sample(nrow(fit_oi_all), 5000)

plot1 <- ggplot(fit_oi_all, aes(x = group, y = beta, group = group, color =group)) +
  geom_quasirandom() +
  geom_boxplot(width = 0.2, color = "black", outliers = F)+
  xlab("")+
  ggtitle("Age-associated TEs")+#, " (", nrow(for_ggplot), " hits)"))+
  ylab("Effect size of age:breed size") +
  theme_blaise +
  # xlim(-1,1)+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.position = 'none'
  )+
  scale_color_manual(
    values = c("red", "darkorange", "lightgreen", "lightblue", "pink", "black"),
    limits = c("Small", "Medium", "Standard", "Large", "Giant")
  )+
  scale_x_discrete(expand = c(0.2,0.01), limits = c("Small", "Medium", "Standard", "Large", "Giant"))+
  scale_y_continuous(expand = c(0.01,0.01))+
  coord_cartesian(clip = "off")+
  geom_signif(
    comparisons = list(c("Small", "Medium"), c("Small", "Standard"), c("Small", "Large"), c("Small", "Giant"),
                       c("Medium", "Large"), c("Medium", "Giant"),
                       c("Standard", "Large"), c("Standard", "Giant"), 
                       c("Large", "Giant")
                       ),
    map_signif_level = T, textsize = 4, na.rm = T, test = 't.test', color = 'black', tip_length = 0, #step_increase = 0.05, 
    y_position = c(0.02,0.0275,0.035,0.0425,0.05,0.0575,0.065,0.0725,0.08,0.0875), vjust = 0.7
  )+
  geom_signif(
    comparisons = list(c("Medium", "Standard")),
    map_signif_level = T, textsize = 4, na.rm = T, test = 't.test', color = 'black', tip_length = 0, #step_increase = 0.05, 
      y_position = c(0.08)
  )+
  labs(tag = 'A')

plot1

t.test(fit_oi_all$beta[fit_oi_all$group == "Small"], fit_oi_all$beta[fit_oi_all$group == "Giant"])
t.test(fit_oi_all$beta[fit_oi_all$group == "Large"], fit_oi_all$beta[fit_oi_all$group == "Giant"])

mean(fit_oi_all$beta[fit_oi_all$group == "Giant"])
mean(fit_oi_all$beta[fit_oi_all$group == "Small"])

fit <- readRDS("../pqlseq_results/pqlseq_res_sexintage_maxgap250.rds")
fit <- fit[fit$converged_Age_at_sample_int_Female == T,]
fit <- fit[fit$converged_Age_at_sample_int_Male == T,]
pqlseq_res <- readRDS("../pqlseq_results/pqlseq_res_maxgap250.rds")
pqlseq_res <- pqlseq_res[order(pqlseq_res$padj),]
pqlseq_res <- pqlseq_res[pqlseq_res$converged == T,]

pqlseq_res$count <- 1:nrow(pqlseq_res)
pqlseq_res$fdr <- pqlseq_res$count * pqlseq_res$padj
pqlseq_res$fdr_perc <- pqlseq_res$fdr / pqlseq_res$count
fit <- fit[rownames(fit) %in% rownames(pqlseq_res)[pqlseq_res$fdr_perc < 0.01 & pqlseq_res$beta < 0],]

col_oi <- "TE"
for (col in c("Female","Male")){
  
  fit_ <- fit[,grepl(col, colnames(fit))]
  
  new_colname <-  paste0("Age_at_sample_int_", gsub(".*Sex","",col))
  
  colnames(fit_) <- gsub(paste0("_",new_colname), "", colnames(fit_))
  
  fit_ <- fit_[order(fit_$padj),]
  fit_ <- fit_[fit_$converged == T,]
  
  fit_$count <- 1:nrow(fit_)
  fit_$fdr <- fit_$count * fit_$padj
  fit_$fdr_perc <- fit_$fdr / fit_$count
  
  age_effect_plots_fdr05=list()
  region_metaData_oi <- region_metaData[region_metaData[,c(paste("TE"))] == 1,]
  fit_oi <- fit_[rownames(fit_) %in% region_metaData_oi$region,]
  
  fit_oi$group <- col
  
  if (col == "Female") { fit_oi_all <- fit_oi} else {fit_oi_all <- rbind(fit_oi_all, fit_oi)}
}

# fit_oi_all$beta <- abs(fit_oi_all$beta)

plot2 <- ggplot(fit_oi_all, aes(x = group, y = beta, group = group, color =group)) +
  geom_quasirandom() +
  geom_boxplot(width = 0.2, color = "black", outliers = F)+
  # geom_density(data = fit_oi_all[fit_oi_all$group == "Male",],alpha = 0.2, linewidth = 1) +
  # geom_density(data = fit_oi_all[fit_oi_all$group == "Female",],alpha = 0.2, linewidth = 1) +
  # geom_histogram(fill = paste0(density_plot_color), color = "black", alpha = 0.5)+
  ylab("Effect size of age:sex") +
  xlab("")+#, " (", nrow(for_ggplot), " hits)"))+
  ggtitle("Age-associated TEs") +
  theme_blaise +
  scale_color_manual(
    limits = c("Male", "Female"),
    values = c("red", "lightblue")
  )+
  # xlim(-1,1)+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.position = 'none'
  )+ 
  # scale_x_discrete(expand = c(0.01,0.01))+ 
  scale_y_continuous(expand = c(0.01,0.01))+
  coord_cartesian(clip = "off")+
  geom_signif(
    comparisons = list(c("Male", "Female")),
    map_signif_level = TRUE, textsize = 5, na.rm = T, test = 't.test',  step_increase = 0.1, color = "black", tip_length = 0, #vjust = 0.6
  )+
  labs(tag = 'B')

plot2

###heritability
#is DNAm in age-associated DMRs more heritable than in non-age-associated DMRs?
fit <- readRDS("../pqlseq_results/pqlseq_res_maxgap250.rds")
fit <- fit[fit$converged == T,]
fit <- na.omit(fit)
fit$chr <- gsub("_.*", "", rownames(fit))
fit$chr <- gsub("chr", "", fit$chr)
fit$count <- 1:nrow(fit)
fit$fdr <- fit$count * fit$padj
fit$fdr_perc <- fit$fdr / fit$count
fit_Y <- fit[grepl("chrY",rownames(fit)),]
fit <- fit[!grepl("chrY", rownames(fit)),]

level_order = c(paste0(c(1:38, "X")))

region_metaData_ <- region_metaData[region_metaData$region %in% rownames(fit),]
region_metaData_ <- region_metaData_[order(region_metaData_$region),]
fit <- fit[order(rownames(fit)),]
all(region_metaData_$region == rownames(fit))

fit$group[fit$fdr_perc < 0.05] <- "age-associated"
fit$group[fit$fdr_perc > 0.05] <- "not age-associated"

# fit <- fit[fit$h2>0,]

plot3 <- ggplot(fit, aes(x = factor(group, levels = c("not age-associated", "age-associated")), y = h2)) +
  geom_quasirandom(color = "orange") +
  geom_boxplot(width = 0.2, color = "black", outliers = F)+
  # geom_density(data = fit[fit$TE_label == "TE",], fill = "purple")+
  # geom_density(data = fit[fit$TE_label == "notTE",], fill = "grey")+
  ylab("heritability") + 
  xlab("")+
  theme_blaise+
  labs(tag = 'C') + 
  theme(legend.position = "none")+
  coord_cartesian(clip = "off")+
  geom_signif(
    comparisons = list(c("not age-associated", "age-associated")),
    map_signif_level = TRUE, textsize = 5, na.rm = T, test = 't.test',  step_increase = 0.1, color = "black", tip_length = 0
  )#+
# geom_label_repel(aes(label=ifelse(label!="0",
#                                   as.character(label),'')),
#                  box.padding = .3,
#                  point.padding = 0.01, segment.color = 'darkgrey', color = "black", size = 6
# )

####
fit <- readRDS("../pqlseq_results/pqlseq_res_maxgap250.rds")
fit <- fit[fit$converged == T,]
fit <- na.omit(fit)
fit$chr <- gsub("_.*", "", rownames(fit))
fit$chr <- gsub("chr", "", fit$chr)
fit$count <- 1:nrow(fit)
fit$fdr <- fit$count * fit$padj
fit$fdr_perc <- fit$fdr / fit$count
fit_Y <- fit[grepl("chrY",rownames(fit)),]
fit <- fit[!grepl("chrY", rownames(fit)),]

level_order = c(paste0(c(1:38, "X")))

region_metaData_ <- region_metaData[region_metaData$region %in% rownames(fit),]
region_metaData_ <- region_metaData_[order(region_metaData_$region),]
fit <- fit[order(rownames(fit)),]
all(region_metaData_$region == rownames(fit))

#####genes####
fit$gene_id <- as.character(region_metaData_$gene_id)
fit$Promoter_id <- as.character(region_metaData_$Promoter_id)
fit$Exon <- as.character(region_metaData_$exon)
fit$TE <- as.character(region_metaData_$LINE)

# fit <- fit[fit$h2>0,]
fit <- fit[fit$fdr_perc<0.05,]

fit$group <- "notGenebody"
fit$group[fit$gene_id != '0'] <- 'Genebody'

fit_ <- fit
fit_$group <- "notPromoter"
fit_$group[fit$Promoter_id != '0'] <- 'Promoter'
fit <- rbind(fit, fit_)

fit_ <- fit
fit_$group <- "notTE"
fit_$group[fit$TE != '0'] <- 'TE'
fit <- rbind(fit, fit_)
plot4 <- ggplot(fit, aes(x = factor(group, levels = c("notGenebody", "Genebody", "notPromoter", "Promoter","notTE","TE")), 
                         y = h2, color = group)) +
  geom_quasirandom() +
  geom_boxplot(width = 0.2, color = "black", outliers = F)+
  # geom_density(data = fit[fit$TE_label == "TE",], fill = "purple")+
  # geom_density(data = fit[fit$TE_label == "notTE",], fill = "grey")+
  ylab("heritability") + 
  xlab("")+
  ggtitle("Heritability of age-associated regions")+
  theme_blaise+
  labs(tag = 'D') +
  theme(legend.position = "none")+
  coord_cartesian(clip = "off")+
  geom_signif(
    comparisons = list(c("notGenebody", "Genebody"), c("notPromoter", "Promoter"), c("notTE", "TE")),
    map_signif_level = T, textsize = 5, na.rm = T, test = 't.test',  step_increase = 0, color = "black", tip_length = 0
  )#+
# geom_label_repel(aes(label=ifelse(gene_id!="0" & h2>0.3,
#                                   as.character(gene_id),'')),
#                  box.padding = .3,
#                  point.padding = 0.01, segment.color = 'darkgrey', color = "black", size = 6
# )

plot5 <- plot3 + plot4 + plot_layout(widths = c(1,2))
plot5

svglite("DAP_FIGURE3.svg", fix_text_size = F, width = 14, height = 10)
print((plot1+plot2) / plot5)
dev.off()
