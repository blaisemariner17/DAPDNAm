#'plot the odds ratio scatter plot for enriched and depleted regions from the hyper and hypo methylated regions of interest
#'
#' @param pqlseq_res pqlseq_results
#' @param region_metaData region metadata
#' @param fdr_perc fdr percent cutoff
#' @param omit_class what classes to omit from the scatterplot. e.g. c("SINE","LINE")
#' @param label_order the labels to be assigned to the plots in ggplot
#' @return Function returns dataframe of oddsratio calculations for hypo and hyper comparisons of interest
#' @export odds_ratio_hyper_hypo_plot

odds_ratio_hyper_hypo_plot <- function(pqlseq_res,
                                       region_metaData = region_metaData,
                                       fdr_perc = 0.05,
                                       omit_class = c(),
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

) {
  pqlseq_res <- pqlseq_res[pqlseq_res$converged == T,]
  region_metaData <- region_metaData[region_metaData$region %in% rownames(pqlseq_res),]

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
  region_metaData_fdr1_oi <- region_metaData_fdr1[region_metaData_fdr1$region %in% rownames(pqlseq_res_fdr1)[pqlseq_res_fdr1$beta > 0],]
  for (col in colnames(region_metaData_fdr1_oi)[!colnames(region_metaData_fdr1_oi) %in% omit_class]){
    if (col == "distance_nearest_gene_start") {
      next
    }
    if (1 %in% region_metaData_fdr1_oi[,c(paste(col))] &
        0 %in% region_metaData_fdr1_oi[,c(paste(col))]) {
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
      print(data_oi)
      odds_ratio <- epitools::oddsratio(data_oi)
      #odds_ratio

      # if (odds_ratio$p.value[2,2] < 0.05) { color = "red"} else {color = "black"}

      res_ <- data.frame("class" = col, "pval" = round(odds_ratio$p.value[2,2], digits = 5), "odds_ratio_log2" = round(log2(odds_ratio$measure[2,1]), digits = 5),
                         "lower" = round(log2(odds_ratio$measure[2,2]), digits = 5), "upper" = round(log2(odds_ratio$measure[2,3]), digits = 5))
      if (col == col in colnames(region_metaData_fdr1_oi)[!colnames(region_metaData_fdr1_oi) %in% omit_class][1]){
        for_ggplot_hyper <- res_
      } else {
        for_ggplot_hyper <- rbind(for_ggplot_hyper, res_)
      }
    } else {next}
  }

  #hypo
  region_metaData_fdr1_oi <- region_metaData_fdr1[region_metaData_fdr1$region %in% rownames(pqlseq_res_fdr1)[pqlseq_res_fdr1$beta < 0],]
  for (col in colnames(region_metaData_fdr1_oi)[!colnames(region_metaData_fdr1_oi) %in% omit_class]){
    if (col == "distance_nearest_gene_start") {
      next
    }
    if (1 %in% region_metaData_fdr1_oi[,c(paste(col))] &
        0 %in% region_metaData_fdr1_oi[,c(paste(col))]) {
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
      print(data_oi)
      odds_ratio <- epitools::oddsratio(data_oi)
      #odds_ratio

      res_ <- data.frame("class" = col, "pval" = round(odds_ratio$p.value[2,2], digits = 2), "odds_ratio_log2" = round(log2(odds_ratio$measure[2,1]), digits = 5),
                         "lower" = round(log2(odds_ratio$measure[2,2]), digits = 5), "upper" = round(log2(odds_ratio$measure[2,3]), digits = 5))
      if (col == colnames(region_metaData_fdr1_oi)[!colnames(region_metaData_fdr1_oi) %in% omit_class][1]){
        for_ggplot_hypo <- res_
      } else {
        for_ggplot_hypo <- rbind(for_ggplot_hypo, res_)
      }
    } else {next}
  }

  #padj
  for_ggplot_hypo$padj_hypo <- round(p.adjust(for_ggplot_hypo$pval,method="BH"),digits = 10)
  for_ggplot_hyper$padj_hyper <- round(p.adjust(for_ggplot_hyper$pval,method="BH"),digits = 10)

  #oddsratio log2
  for_ggplot_hyper$odds_ratio_log2_hyper <- for_ggplot_hyper$odds_ratio_log2
  for_ggplot_hypo$odds_ratio_log2_hypo <- for_ggplot_hypo$odds_ratio_log2
  #class
  for_ggplot_hypo$class_hypo <- for_ggplot_hypo$class
  for_ggplot_hyper$class_hyper <- for_ggplot_hyper$class
  #and rownames assignment
  rownames(for_ggplot_hypo) <- for_ggplot_hypo$class
  rownames(for_ggplot_hyper) <- for_ggplot_hyper$class
  #get the classes that are in the other
  for_ggplot_hypo <- for_ggplot_hypo[rownames(for_ggplot_hypo) %in% rownames(for_ggplot_hyper),]
  for_ggplot_hyper <- for_ggplot_hyper[rownames(for_ggplot_hyper) %in% rownames(for_ggplot_hypo),]

  for_ggplot_hyper <- for_ggplot_hyper[order(for_ggplot_hyper$class),]
  for_ggplot_hypo <- for_ggplot_hypo[order(for_ggplot_hypo$class),]

  if(all(rownames(for_ggplot_hypo) == rownames(for_ggplot_hyper))==F){return(message("hypo results do not match hyper results"))}

  odds_ratio_results <- cbind(for_ggplot_hyper[,c("class_hyper", "padj_hyper", "odds_ratio_log2_hyper")],
                           for_ggplot_hypo[,c("class_hypo", "padj_hypo", "odds_ratio_log2_hypo")])

  odds_ratio_results$significance_label <- "No significance (p<sub>adj</sub> > 0.05, Fisher Exact test)"

  odds_ratio_results$significance_label[odds_ratio_results$padj_hyper < 0.05 & odds_ratio_results$odds_ratio_log2_hyper > 0] <- "Enriched in hyper only"
  odds_ratio_results$significance_label[odds_ratio_results$padj_hypo < 0.05 & odds_ratio_results$odds_ratio_log2_hypo < 0] <- "Depleted in hypo only"
  odds_ratio_results$significance_label[odds_ratio_results$padj_hyper < 0.05 & odds_ratio_results$odds_ratio_log2_hyper < 0] <- "Depleted in hyper only"
  odds_ratio_results$significance_label[odds_ratio_results$padj_hypo < 0.05 & odds_ratio_results$odds_ratio_log2_hypo > 0] <- "Enriched in hypo only"
  odds_ratio_results$significance_label[odds_ratio_results$padj_hyper < 0.05 & odds_ratio_results$padj_hypo < 0.05 &
                                       (odds_ratio_results$odds_ratio_log2_hyper < 0 & odds_ratio_results$odds_ratio_log2_hypo > 0)] <- "Enriched in hypo & depleted in hyper"
  odds_ratio_results$significance_label[odds_ratio_results$padj_hyper < 0.05 & odds_ratio_results$padj_hypo < 0.05 &
                                       (odds_ratio_results$odds_ratio_log2_hyper >0 & odds_ratio_results$odds_ratio_log2_hypo < 0)] <- "Enriched in hyper & depleted in hypo"
  odds_ratio_results$significance_label[odds_ratio_results$padj_hyper < 0.05 & odds_ratio_results$padj_hypo < 0.05 &
                                       odds_ratio_results$odds_ratio_log2_hyper >0 & odds_ratio_results$odds_ratio_log2_hypo >0] <- "Enriched in hyper & hypo"
  odds_ratio_results$significance_label[odds_ratio_results$padj_hyper < 0.05 & odds_ratio_results$padj_hypo < 0.05 &
                                       odds_ratio_results$odds_ratio_log2_hyper <0 & odds_ratio_results$odds_ratio_log2_hypo <0] <- "Depleted in hyper & hypo"

  odds_ratio_results$significance_label <- factor(odds_ratio_results$significance_label,
                                               levels = label_order)

  odds_ratio_results <- odds_ratio_results[order(odds_ratio_results$significance_label),]

  return(odds_ratio_results)
}
