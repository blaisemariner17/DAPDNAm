#'plot the odds ratio scatter plot for enriched and depleted regions from the hyper and hypo methylated regions of interest
#'
#' @param pqlseq_res pqlseq_results
#' @param region_metaData region metadata
#' @param fdr_perc fdr percent cutoff
#' @param omit_class what classes to omit from the scatterplot. e.g. c("SINE","LINE")
#' @param color_order the colors to be assigned to the plots in ggplot
#' @return Function returns dataframe of oddsratio calculations for hypo and hyper comparisons of interest
#' @export odds_ratio_hyper_hypo_plot

odds_ratio_hyper_hypo_plot <- function(pqlseq_res,
                                       region_metaData = region_metaData,
                                       fdr_perc = 0.05,
                                       omit_class = c()
) {
  pqlseq_res <- pqlseq_res[pqlseq_res$converged == T,]
  pqlseq_res <- pqlseq_res[order(pqlseq_res$padj),]
  pqlseq_res$name <- rownames(pqlseq_res)

  pqlseq_res$count <- 1:nrow(pqlseq_res)
  pqlseq_res$fdr <- pqlseq_res$count * pqlseq_res$padj
  pqlseq_res$fdr_perc <- pqlseq_res$fdr / pqlseq_res$count
  pqlseq_res_fdr1 <- pqlseq_res[pqlseq_res$fdr_perc < fdr_perc,]
  region_metaData_fdr1 <- region_metaData[region_metaData$region %in% rownames(pqlseq_res_fdr1),]
  pqlseq_res_fdr1 <- pqlseq_res_fdr1[order(rownames(pqlseq_res_fdr1)),]
  region_metaData_fdr1 <- region_metaData_fdr1[order(region_metaData_fdr1$region),]

  if (all(region_metaData_fdr1$region == rownames(pqlseq_res_fdr1)) == FALSE) {return(message("region_metadata does not match pqlseq_res"))}

  #pos beta
  if (exists("for_ggplot_hyper")){rm(for_ggplot_hyper)}
  res <- list()
  region_metaData_fdr1_oi <- region_metaData_fdr1[region_metaData_fdr1$region %in% rownames(pqlseq_res_fdr1)[pqlseq_res_fdr1$beta > 0],]
  for (col in colnames(region_metaData_fdr1_oi)[!colnames(region_metaData_fdr1_oi) %in% omit_class]){
    if (col == "distance_nearest_gene_start") {
      next
    }
    if (1 %in% region_metaData_fdr1_oi[,c(paste(col))]){
      class_oi <- c(paste0(col), paste0("not_", col))
      outcome_oi <- c("Hypermethylated with Age", "Not")

      number_sighyp_from_col <- nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 1,])
      number_notsighyp_from_col <- sum(region_metaData[,c(paste0(col))] == 1) - nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 1,])
      number_sig_notfrom_col <- nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 0,])
      numer_notsig_notfrom_col <- sum(region_metaData[,c(paste0(col))] == 0) - nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 0,])

      data_oi <- matrix(c(number_sighyp_from_col,
                          number_sig_notfrom_col,
                          number_notsighyp_from_col,
                          numer_notsig_notfrom_col),
                        nrow=2, ncol=2, byrow=TRUE)

      colnames(data_oi) <- class_oi
      rownames(data_oi) <- outcome_oi
      data_oi
      odds_ratio <- oddsratio(data_oi)
      odds_ratio

      # if (odds_ratio$p.value[2,2] < 0.05) { color = "red"} else {color = "black"}

      res_ <- data.frame("class" = col, "pval" = round(odds_ratio$p.value[2,2], digits = 2), "odds_ratio_log10" = round(log10(odds_ratio$measure[2,1]), digits = 3),
                         "lower" = round(log10(odds_ratio$measure[2,2]), digits = 3), "upper" = round(log10(odds_ratio$measure[2,3]), digits = 3))
      res[[col]] <- res_
      if (exists("for_ggplot_hyper") == FALSE){
        for_ggplot_hyper <- res_
      } else {
        for_ggplot_hyper <- rbind(for_ggplot_hyper, res_)
      }
    }
  }
  for_ggplot_hyper$padj <- round(p.adjust(for_ggplot_hyper$pval,method="BH"),digits = 4)
  for_ggplot_hyper$color <- "black"
  for_ggplot_hyper$color[for_ggplot_hyper$padj < 0.05] <- "red"
  for_ggplot_hyper$class <- factor(for_ggplot_hyper$class,                                    # Factor levels in decreasing order
                                   levels = for_ggplot_hyper$class[order(for_ggplot_hyper$odds_ratio_log10, decreasing = FALSE)])

  #hypo
  if (exists("for_ggplot_hypo")){rm(for_ggplot_hypo)}
  res <- list()
  region_metaData_fdr1_oi <- region_metaData_fdr1[region_metaData_fdr1$region %in% rownames(pqlseq_res_fdr1)[pqlseq_res_fdr1$beta < 0],]
  for (col in colnames(region_metaData_fdr1_oi)colnames(region_metaData_fdr1_oi)[!colnames(region_metaData_fdr1_oi) %in% omit_class]){
    if (col == "distance_nearest_gene_start") {
      next
    }
    if (1 %in% region_metaData_fdr1_oi[,c(paste(col))]){
      class_oi <- c(paste0(col), paste0("not_", col))
      outcome_oi <- c("Hypomethylated with Age", "Not")

      number_sighyp_from_col <- nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 1,])
      number_notsighyp_from_col <- sum(region_metaData[,c(paste0(col))] == 1) - nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 1,])
      number_sig_notfrom_col <- nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 0,])
      numer_notsig_notfrom_col <- sum(region_metaData[,c(paste0(col))] == 0) - nrow(region_metaData_fdr1_oi[region_metaData_fdr1_oi[,c(paste0(col))] == 0,])

      data_oi <- matrix(c(number_sighyp_from_col,
                          number_sig_notfrom_col,
                          number_notsighyp_from_col,
                          numer_notsig_notfrom_col),
                        nrow=2, ncol=2, byrow=TRUE)

      colnames(data_oi) <- class_oi
      rownames(data_oi) <- outcome_oi
      data_oi
      odds_ratio <- oddsratio(data_oi)
      odds_ratio

      res_ <- data.frame("class" = col, "pval" = round(odds_ratio$p.value[2,2], digits = 2), "odds_ratio_log10" = round(log10(odds_ratio$measure[2,1]), digits = 3),
                         "lower" = round(log10(odds_ratio$measure[2,2]), digits = 3), "upper" = round(log10(odds_ratio$measure[2,3]), digits = 3))
      res[[col]] <- res_
      if (exists("for_ggplot_hypo") == FALSE){
        for_ggplot_hypo <- res_
      } else {
        for_ggplot_hypo <- rbind(for_ggplot_hypo, res_)
      }
    }
  }

  for_ggplot_hypo$padj <- round(p.adjust(for_ggplot_hypo$pval,method="BH"),digits = 4)
  for_ggplot_hypo$color <- "black"
  for_ggplot_hypo$color[for_ggplot_hypo$padj < 0.05] <- "red"
  for_ggplot_hypo$class <- factor(for_ggplot_hypo$class,                                    # Factor levels in decreasing order
                                  levels = for_ggplot_hypo$class[order(for_ggplot_hypo$odds_ratio_log10, decreasing = FALSE)])

  for_ggplot_hyper$odds_ratio_log10_hyper <- for_ggplot_hyper$odds_ratio_log10
  for_ggplot_hyper$class_hyper <- for_ggplot_hyper$class
  for_ggplot_hyper$color_hyper <- for_ggplot_hyper$color
  rownames(for_ggplot_hyper) <- for_ggplot_hyper$class

  for_ggplot_hypo$odds_ratio_log10_hypo <- for_ggplot_hypo$odds_ratio_log10
  for_ggplot_hypo$class_hypo <- for_ggplot_hypo$class
  for_ggplot_hypo$color_hypo <- for_ggplot_hypo$color
  rownames(for_ggplot_hypo) <- for_ggplot_hypo$class

  for_ggplot_hypo <- for_ggplot_hypo[rownames(for_ggplot_hypo) %in% rownames(for_ggplot_hyper),]
  for_ggplot_hyper <- for_ggplot_hyper[rownames(for_ggplot_hyper) %in% rownames(for_ggplot_hypo),]

  if(all(rownames(for_ggplot_hypo) == rownames(for_ggplot_hyper))==F){return(message("hypo results do not match hyper results"))}

  for_ggplot_both <- cbind(for_ggplot_hyper[,c("class_hyper", "color_hyper", "odds_ratio_log10_hyper")],
                           for_ggplot_hypo[,c("class_hypo", "color_hypo", "odds_ratio_log10_hypo")])

  for_ggplot_both$color <- "black"

  for_ggplot_both$color[for_ggplot_both$color_hyper == "red" & for_ggplot_both$odds_ratio_log10_hyper > 0] <- "blue"
  for_ggplot_both$color[for_ggplot_both$color_hypo == "red" & for_ggplot_both$odds_ratio_log10_hypo < 0] <- "darkgreen"
  for_ggplot_both$color[for_ggplot_both$color_hyper == "red" & for_ggplot_both$odds_ratio_log10_hyper < 0] <- "darkorange"
  for_ggplot_both$color[for_ggplot_both$color_hypo == "red" & for_ggplot_both$odds_ratio_log10_hypo > 0] <- "maroon"
  for_ggplot_both$color[for_ggplot_both$color_hyper == "red" & for_ggplot_both$color_hypo == "red" &
                          (for_ggplot_both$odds_ratio_log10_hyper < 0 & for_ggplot_both$odds_ratio_log10_hypo > 0)] <- "red"
  for_ggplot_both$color[for_ggplot_both$color_hyper == "red" & for_ggplot_both$color_hypo == "red" &
                          (for_ggplot_both$odds_ratio_log10_hyper >0 & for_ggplot_both$odds_ratio_log10_hypo < 0)] <- "purple"

  for_ggplot_both$color <- factor(for_ggplot_both$color,
                                  levels =color_order)

  for_ggplot_both <- for_ggplot_both[order(for_ggplot_both$color),]

  # theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 24, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"), axis.text = element_text(family = "sans", size = 18, color="black"),axis.title.y = element_markdown(family = "sans", size = 20), axis.title.x = element_markdown(family = "sans", size = 20),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()
  #
  # plot_ <- ggplot(for_ggplot_both[!for_ggplot_both$class_hyper %in% omit_class,],
  #                 aes(x = odds_ratio_log10_hypo, y = odds_ratio_log10_hyper, label = class_hyper, color = color)) +
  #   geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") +
  #   geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") +
  #   geom_point(size = 3) +
  #   ylab("hypermethylation log<sub>10</sub>(odds ratio)")+
  #   xlab("hypomethylation log<sub>10</sub>(odds ratio)")+
  #   theme_blaise +
  #   theme(plot.caption = element_markdown(size = 14), legend.position = c(0.8, 0.82),
  #         legend.key.height = unit(1, 'line'),
  #         legend.key.width = unit(1, 'line'),
  #         legend.text = element_markdown(size = 12)
  #   )+
  #   geom_label(x = -0.25, y = 0.65, label = "Significantly enriched in hypermethylated regions \n and significantly depleted in hypomethylated regions",
  #              family= 'sans', color = "purple") +
  #   geom_label(x = 0.4, y = -0.65, label = "Significantly enriched in hypomethylated regions \n and significantly depleted in hypermethylated regions",
  #              family= 'sans', color = "red") +
  #   geom_label_repel(aes(label=ifelse(color!="black", as.character(class_hyper),'')), box.padding = .2,
  #                    point.padding = 0.01, segment.color = 'grey', color = "black"
  #   ) +
  #   xlim(-0.5,0.6) + ylim(-0.82,0.77) +
  #   scale_color_manual(values = color_order,
  #                      labels = c(
  #                        "Enriched in hyper & depleted in hypo",
  #                        "Enriched in hyper only",
  #                        "Depleted in hypo only",
  #                        "No significance (p<sub>adj</sub> > 0.05, Fisher Exact test)",
  #                        "Depleted in hyper only",
  #                        "Enriched in hypo only",
  #                        "Enriched in hypo & depleted in hyper"
  #                      ))

  return(for_ggplot_both)
}
