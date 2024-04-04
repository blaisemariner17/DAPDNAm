#### Blaise Mariner now has a slight grasp on Sol and the dap mDNA data so now it is time to make a full workflow pipeline
#### this is part 2 of this series where we can visualize the regions, their coverage and methylation that we made from pt 1. 
#### part 3 will filter out the samples with a low lib size
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-01-12
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

col_oi = c(
  "Promoter", "gene_bool","exon", "intron",
  "CpG_island", "CpG_shore", "CpG_shelf",
  # "ChrSt_promoter", 
  "ChrSt_quies","ChrSt_heterochromatin",
  # "ChrSt_polycomb", 
  "ChrSt_enhancer",
  # "DNA.transposon", "Retrotransposon",
  "TE"
)
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
  "patchwork"
)
lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 12, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      
                      axis.text = element_text(family = "sans", size = 14, color="black"),#axis.title.y = element_markdown(family = "sans", size = 14),   
                      axis.title.x = element_markdown(family = "sans", size = 14),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

#1.5 hrs ish with 25 cores (1 managing) and 250 Gb mem (200 errors out) during peak hours
cores_ <- 24

# this sets the working directory to this script's path
# if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
# }
print(getwd())

# metadata
metaData <- read_rds("../metadata_samples/P1-DAP-metaData-2400216.rds")
p1 <- metaData

### let's look at some effect sizes ####
fit <- read_rds("../pqlseq_results/pqlseq_res_maxgap250.rds")
summary(fit)
fit <- fit[fit$converged == T,]
fit <- fit[order(fit$padj),]
#####
age_effect_plots=list()
for (col in col_oi){
  region_metaData_oi <- region_metaData[region_metaData[,c(paste(col))] == 1,]
  fit_oi <- fit[rownames(fit) %in% region_metaData_oi$region,]
  fit_oi$beta <- fit_oi$beta / max(fit_oi$beta)
  
  for_ggplot <- data.frame(fit_oi$beta)
  colnames(for_ggplot) <- c("data")
  ylab_ = col
  if (col == "Promoter"){ylab_ = "promoter"}
  if (col == "gene_bool"){ylab_ = "gene body"}
  if (col == "CpG_island"){ylab_ = "CpG island"}
  if (col == "CpG_shore"){ylab_ = "CpG shore"}
  if (col == "CpG_shelf"){ylab_ = "CpG shelf"}
  if (col == "ChrSt_enhancer"){ylab_ = "enhancer"}
  if (col == "ChrSt_polycomb"){ylab_ = "polycomb"}
  if (col == "ChrSt_quies"){ylab_ = "quiescent"}
  if (col == "ChrSt_heterochromatin"){ylab_ = "heterochromatin"}
  plot_ <- ggplot(for_ggplot, aes(x = data)) +
    geom_density(fill = "lightgrey") +
    xlab("")+
    ylab(paste0(ylab_))+#, " (n = ", nrow(for_ggplot), ")"))+
    ggtitle("") +
    theme_blaise +
    xlim(-0.45,0.45)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), plot.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          axis.title.y = element_text(angle = 0,size=18, hjust = 1, vjust = 0.5)            
    )+coord_cartesian(clip = "off")+
    geom_vline(linetype = "dotted", color = "red", xintercept = 0, linewidth = 1)
  # plot_
  if (col == col_oi[length(col_oi)]) {
    plot_ <- plot_ + xlab("Effect size of age") + theme(axis.text.x=element_text(family = "sans", size = 18, color="black"),
                                                        axis.line.x = element_line(colour = "black", linewidth = 1),
                                                        axis.title.x = element_text(family = "sans", size = 18, hjust = 0.5, color="black")
    )
  } else {
    plot_ <- plot_ + theme(axis.line.x = element_blank())
  }
  
  age_effect_plots[[col]] <- plot_
}
age_effect_plots[["Promoter"]] <- age_effect_plots[["Promoter"]]  + labs(tag = 'A')

plot1 <- (patchwork::wrap_plots(age_effect_plots, nrow = length(age_effect_plots), ncol = 1))

fit$count <- 1:nrow(fit)
fit$fdr <- fit$count * fit$padj
fit$fdr_perc <- fit$fdr / fit$count

age_effect_plots_fdr05=list()
for (col in col_oi){
  region_metaData_oi <- region_metaData[region_metaData[,c(paste(col))] == 1,]
  fit_oi <- fit[rownames(fit) %in% region_metaData_oi$region & fit$fdr_perc < 0.05,]
  fit_oi$beta <- fit_oi$beta / max(fit_oi$beta)
  
  for_ggplot <- data.frame(fit_oi$beta)
  if(nrow(for_ggplot) == 0){next}
  colnames(for_ggplot) <- c("data")
  
  plot_ <- ggplot(for_ggplot, aes(x = data)) +
    geom_density(fill = "lightgrey") +
    xlab("")+
    ylab(paste0(col, " (", nrow(for_ggplot), " hits)"))+
    ggtitle("") +
    theme_blaise +
    xlim(-0.45,0.45)+
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), plot.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          axis.title.y = element_blank(),
          axis.line.y=element_blank()
    )+
    geom_vline(linetype = "dotted", color = "red", xintercept = 0, linewidth = 1)
  
  if (col == col_oi[length(col_oi)]) {
    plot_ <- plot_ + xlab("Effect size of age") + theme(axis.text.x=element_text(family = "sans", size = 18, color="black"),
                                                        axis.line.x = element_line(colour = "black", linewidth = 1),
                                                        axis.title.x = element_text(family = "sans", size = 18, hjust = 0.5, color="black")
    )
  } else {
    plot_ <- plot_ + theme(axis.line.x = element_blank())
  }
  
  age_effect_plots_fdr05[[col]] <- plot_
}

age_effect_plots_fdr05[["Promoter"]] <- age_effect_plots_fdr05[["Promoter"]] # + labs(tag = 'B')

plot2 <- (patchwork::wrap_plots(age_effect_plots_fdr05, nrow = length(age_effect_plots_fdr05), ncol = 1))

fit <- readRDS("../pqlseq_results/pqlseq_res_maxgap250.rds")
fit <- fit[fit$converged == T,]
fit$count <- 1:nrow(fit)
fit$fdr <- fit$count * fit$padj
fit$fdr_perc <- fit$fdr / fit$count
age_effect_plots_fdr05_bar=list()
for (col in col_oi){
  region_metaData_oi <- region_metaData[region_metaData[,c(paste(col))] == 1,]
  fit_oi <- fit[rownames(fit) %in% region_metaData_oi$region & fit$fdr_perc < 0.05,]

  
  total <- nrow(fit[rownames(fit) %in% region_metaData_oi$region,])
  fdr05_hyper <- nrow(fit_oi[fit_oi$fdr_perc < 0.05 & fit_oi$beta > 0,])
  fdr05_hypo <- nrow(fit_oi[fit_oi$fdr_perc < 0.05 & fit_oi$beta < 0,])
  
  for_ggplot <- data.frame("ann" = c("Total regions", 
                                     "Hypermethylated regions",
                                     "Hypomethylated regions"), 
                           "ann_short" = c("Total", 
                                           "Hypermethylated",
                                           "Hypomethylated"),
                           "Number" = c(total, fdr05_hyper, fdr05_hypo))
  for_ggplot$ann_short_total <- ''
  for_ggplot$ann_short_total[for_ggplot$ann_short == "Total"] <- paste0("",#for_ggplot$ann_short[for_ggplot$ann_short == "Total"], 
                                                                        "", for_ggplot$Number[for_ggplot$ann_short == "Total"])
  for_ggplot$ann_short_hyp <- ''
  for_ggplot$ann_short_hyp[for_ggplot$ann_short != "Total"] <- paste0(for_ggplot$ann_short[for_ggplot$ann_short != "Total"], ": ", for_ggplot$Number[for_ggplot$ann_short != "Total"])
  plot_ <- ggplot(for_ggplot, aes(x = Number, 
                                  y = factor(ann, levels = c("Hypomethylated regions","Hypermethylated regions", "Total regions")), 
                                  fill = ann)) +
    geom_bar(stat = 'identity')+
    geom_text(aes(label = ann_short_hyp), colour = "black", hjust = -0.01, size = 4)+
    geom_text(aes(label = ann_short_total), colour = "white", hjust = 1, size = 4)+
    xlab("")+
    ylab(paste0(col))+#, " (", nrow(for_ggplot), " hits)"))+
    ggtitle("") +
    theme_blaise +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), plot.title = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          axis.title.y = element_blank(),
          axis.line.y=element_blank(), 
          axis.line.x = element_blank(),
          legend.position = 'none'
    )+
    scale_fill_manual(values = c( "purple","red","black" ))+
    scale_x_continuous(expand = c(0,0), limits = c(0,total+2000))
  
  if (col == col_oi[length(col_oi)]) {
    plot_ <- plot_ + xlab("Number of regions") + theme(axis.text.x=element_blank(),
                                                       axis.line.x = element_line(colour = "black", linewidth = 1),
                                                       axis.title.x = element_text(family = "sans", size = 18, hjust = 0.5, color="black")
    )
  } else {
    plot_ <- plot_ + theme(axis.line.x = element_blank())
  }
  
  age_effect_plots_fdr05_bar[[col]] <- plot_
}
age_effect_plots_fdr05_bar[["Promoter"]] <- age_effect_plots_fdr05_bar[["Promoter"]] # + labs(tag = 'C')

plot3 <- (patchwork::wrap_plots(age_effect_plots_fdr05_bar, nrow = length(age_effect_plots_fdr05_bar), ncol = 1))

fit <- readRDS("../pqlseq_results/pqlseq_res_maxgap250.rds")
fit <- fit[! grepl("chrY", rownames(fit)),]
region_metaData <- readRDS("../metadata_regions/metaData_regions_maxgap_250.rds")

for_ggplot_both <- DAPDNAm::odds_ratio_hyper_hypo_plot(fit, 
                                                       region_metaData = region_metaData,
                                                       omit_class = c())

for_ggplot_both <- for_ggplot_both[for_ggplot_both$class_hyper %in% col_oi,]
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "Promoter"] <- "promoter"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "gene_bool"] <- "gene body"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "CpG_island"] <- "CpG island"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "CpG_shore"] <- "CpG shore"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "CpG_shelf"] <- "CpG shelf"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "ChrSt_enhancer"] <- "enhancer"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "ChrSt_polycomb"] <- "polycomb"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "ChrSt_quies"] <- "quiescent"
for_ggplot_both$class_hyper[for_ggplot_both$class_hyper == "ChrSt_heterochromatin"] <- "heterochromatin"

levels_order <- rev(c("promoter", "exon", "intron", "CpG island", "CpG shore", "CpG shelf", "quiescent", "heterochromatin", 
                      # "polycomb", 
                      "enhancer", "DNA.transposon", "Retrotransposon", "TE"))

# for_ggplot_both$class_hyper <- factor(for_ggplot_both$class_hyper,                                    # Factor levels in decreasing order
#                                              levels = for_ggplot_both$class_hyper[order(for_ggplot_both$odds_ratio_log2_hyper, decreasing = FALSE)])
for_ggplot_both$color <- "black"
for_ggplot_both$color[for_ggplot_both$padj_hyper < 0.05] <- "red"

label_order = c(
  "Enriched in hyper & depleted in hypo",
  "Enriched in hyper only",
  "Depleted in hypo only",
  "No significance (p<sub>adj</sub> > 0.05, Fisher Exact test)",
  "Depleted in hyper only",
  "Enriched in hypo only",
  "Enriched in hypo & depleted in hyper",
  "Enriched in hyper & hypo",
  "Depleted in hyper & hypo"
)

interested_class <- c("promoter", "exon", "intron", "TE", "heterochromatin", 
                      # "polycomb", 
                      "enhancer", "CpG island", "CpG shore", "CpG shelf", "quiescent")

theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 24, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"), axis.text = element_text(family = "sans", size = 18, color="black"),axis.title.y = element_markdown(family = "sans", size = 20), axis.title.x = element_markdown(family = "sans", size = 20),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()
color_order =c("purple", "blue", "darkgreen", "black", "darkorange", "maroon", "red","#66b2b2","lightpink" )

color_order_<-color_order[label_order %in% for_ggplot_both$significance_label[for_ggplot_both$class_hyper %in% interested_class]]

plot4 <- ggplot(for_ggplot_both[for_ggplot_both$class_hyper %in% interested_class,],
                aes(y = odds_ratio_log2_hypo, x = odds_ratio_log2_hyper, label = class_hyper, color = significance_label)) +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") +
  geom_point(size = 4) +
  geom_errorbar(aes(xmax = upper_hyper, xmin = lower_hyper))+
  geom_errorbar(aes(ymin = lower_hypo, ymax = upper_hypo))+
  # ggtitle(paste0("Age-related region enrichment; maxgap=", maxgap))+
  xlab("age-related hypermethylation log<sub>2</sub>(OR)")+
  ylab("age-related hypomethylation log<sub>2</sub>(OR)")+
  theme_blaise +
  theme(plot.caption = element_markdown(size = 14), #legend.position = c(0.8, 0.82),
        legend.key.height = unit(1.2, 'line'),
        legend.key.width = unit(1.2, 'line'),
        legend.text = element_markdown(size = 19),
        # legend.key.size = unit(1.5, 'cm'),
        legend.position = c(0.75, 0.82),
        legend.background = element_rect(fill = "white"),
  )+
  geom_label_repel(aes(label=ifelse(significance_label!="!No significance (p<sub>adj</sub> > 0.05, Fisher Exact test)",
                                    as.character(class_hyper),'')),
                   box.padding = .3,
                   point.padding = 0.01, segment.color = 'darkgrey', color = "black", size = 6
  ) +
  ylim(-1.5,.75) +
  xlim(-1,2.5) +
  scale_color_manual(values = color_order_) + labs(tag = 'B')

#plotly::ggplotly(plot_)
# print(plot4)

# plot4 <- ggplot(for_ggplot_both, aes(y = factor(class_hyper, levels = levels_order), x = odds_ratio_log2_hyper)) + 
#   geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") +
#   geom_errorbar(aes(xmin=lower_hyper, xmax=upper_hyper), width = 0.25)+
#   geom_point(size = 2, color = for_ggplot_both$color) +
#   ylab("") + xlab("Hypermethylated region enrichment log<sub>2</sub>OR")+
#   theme_blaise + labs(tag = 'B')
#   
# for_ggplot_both$color <- "black"
# for_ggplot_both$color[for_ggplot_both$padj_hypo < 0.05] <- "red"
# # for_ggplot_both$class_hyper <- factor(for_ggplot_both$class_hyper,                                    # Factor levels in decreasing order
# #                                       levels = for_ggplot_both$class_hyper[order(for_ggplot_both$odds_ratio_log2_hypo, decreasing = FALSE)])
# 
# plot5 <- ggplot(for_ggplot_both, aes(y = factor(class_hyper, levels = levels_order), x = odds_ratio_log2_hypo)) + 
#   geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") +
#   geom_errorbar(aes(xmin=lower_hypo, xmax=upper_hypo), width = 0.25)+
#   geom_point(size = 2, color = for_ggplot_both$color) +
#   ylab("") + xlab("Hypomethylated region enrichment log<sub>2</sub>OR")+
#   theme_blaise + labs(tag = 'C')

# (plot1 | plot2 | plot3) / (plot4 | plot5)

svglite("DAP_FIGURE2.svg", fix_text_size = F, width = 16, height = 16)
print((plot1 | plot2 | plot3) / (plot4))
dev.off()

