#'plot the odds ratio scatter plot for enriched and depleted regions from the hyper and hypo methylated regions of interest
#'
#' @param pqlseq_res pqlseq_results
#' @param region_metaData region metadata
#' @param fdr_perc fdr percent cutoff
#' @param omit_class what classes to omit from the scatterplot. e.g. c("SINE","LINE")
#' @param label_order the labels to be assigned to the plots in ggplot
#' @return Function returns dataframe of oddsratio calculations for hypo and hyper comparisons of interest. It searches through all the region_metadata columns for columns with 0 and 1 in them and does it automatically
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
  # make sure you only have the regions that have converged
  pqlseq_res <- pqlseq_res[pqlseq_res$converged == T,]
  region_metaData <- region_metaData[region_metaData$region %in% rownames(pqlseq_res),]

  pqlseq_res <- pqlseq_res[order(pqlseq_res$padj),]
  pqlseq_res$name <- rownames(pqlseq_res)

  # generate the fdr percentage for odds ratio calculation
  pqlseq_res$count <- 1:nrow(pqlseq_res)
  pqlseq_res$fdr <- pqlseq_res$count * pqlseq_res$padj
  pqlseq_res$fdr_perc <- pqlseq_res$fdr / pqlseq_res$count
  pqlseq_res_fdr1 <- pqlseq_res[pqlseq_res$fdr_perc < fdr_perc,]
  region_metaData_fdr1 <- region_metaData[region_metaData$region %in% rownames(pqlseq_res_fdr1),]
  pqlseq_res_fdr1 <- pqlseq_res_fdr1[order(rownames(pqlseq_res_fdr1)),]
  region_metaData_fdr1 <- region_metaData_fdr1[order(region_metaData_fdr1$region),]

  if (all(region_metaData_fdr1$region == rownames(pqlseq_res_fdr1)) == FALSE) {return(message("region_metadata does not match pqlseq_res"))}

  #firs, we
  region_metaData_fdr1_oi <- region_metaData_fdr1[region_metaData_fdr1$region %in% rownames(pqlseq_res_fdr1)[pqlseq_res_fdr1$beta > 0],]
  i = 1
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
      # print(data_oi)
      odds_ratio <- epitools::oddsratio(data_oi)
      #odds_ratio

      # if (odds_ratio$p.value[2,2] < 0.05) { color = "red"} else {color = "black"}

      res_ <- data.frame("class" = col, "pval" = round(odds_ratio$p.value[2,2], digits = 5), "odds_ratio_log2" = round(log2(odds_ratio$measure[2,1]), digits = 5),
                         "lower" = round(log2(odds_ratio$measure[2,2]), digits = 5), "upper" = round(log2(odds_ratio$measure[2,3]), digits = 5))
      if (i == 1){
        odds_ratio_results_hyper <- res_
        i = 2
      } else {
        odds_ratio_results_hyper <- rbind(odds_ratio_results_hyper, res_)
      }
    } else {next}
  }

  #hypo
  region_metaData_fdr1_oi <- region_metaData_fdr1[region_metaData_fdr1$region %in% rownames(pqlseq_res_fdr1)[pqlseq_res_fdr1$beta < 0],]
  i = 1
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
      # print(data_oi)
      odds_ratio <- epitools::oddsratio(data_oi)
      #odds_ratio

      res_ <- data.frame("class" = col, "pval" = round(odds_ratio$p.value[2,2], digits = 2), "odds_ratio_log2" = round(log2(odds_ratio$measure[2,1]), digits = 5),
                         "lower" = round(log2(odds_ratio$measure[2,2]), digits = 5), "upper" = round(log2(odds_ratio$measure[2,3]), digits = 5))
      if (i == 1){
        odds_ratio_results_hypo <- res_
        i = 2
      } else {
        odds_ratio_results_hypo <- rbind(odds_ratio_results_hypo, res_)
      }
    } else {next}
  }

  #padj
  odds_ratio_results_hypo$padj_hypo <- round(p.adjust(odds_ratio_results_hypo$pval,method="BH"),digits = 10)
  odds_ratio_results_hyper$padj_hyper <- round(p.adjust(odds_ratio_results_hyper$pval,method="BH"),digits = 10)
  #lower
  odds_ratio_results_hyper$lower_hyper <- odds_ratio_results_hyper$lower
  odds_ratio_results_hypo$lower_hypo <- odds_ratio_results_hypo$lower
  #upper
  odds_ratio_results_hyper$upper_hyper <- odds_ratio_results_hyper$upper
  odds_ratio_results_hypo$upper_hypo <- odds_ratio_results_hypo$upper
  #oddsratio log2
  odds_ratio_results_hyper$odds_ratio_log2_hyper <- odds_ratio_results_hyper$odds_ratio_log2
  odds_ratio_results_hypo$odds_ratio_log2_hypo <- odds_ratio_results_hypo$odds_ratio_log2
  #class
  odds_ratio_results_hypo$class_hypo <- odds_ratio_results_hypo$class
  odds_ratio_results_hyper$class_hyper <- odds_ratio_results_hyper$class
  #and rownames assignment
  rownames(odds_ratio_results_hypo) <- odds_ratio_results_hypo$class
  rownames(odds_ratio_results_hyper) <- odds_ratio_results_hyper$class
  #get the classes that are in the other
  odds_ratio_results_hypo <- odds_ratio_results_hypo[rownames(odds_ratio_results_hypo) %in% rownames(odds_ratio_results_hyper),]
  odds_ratio_results_hyper <- odds_ratio_results_hyper[rownames(odds_ratio_results_hyper) %in% rownames(odds_ratio_results_hypo),]

  odds_ratio_results_hyper <- odds_ratio_results_hyper[order(odds_ratio_results_hyper$class),]
  odds_ratio_results_hypo <- odds_ratio_results_hypo[order(odds_ratio_results_hypo$class),]

  if(all(rownames(odds_ratio_results_hypo) == rownames(odds_ratio_results_hyper))==F){return(message("hypo results do not match hyper results"))}

  odds_ratio_results <- cbind(odds_ratio_results_hyper[,c("class_hyper", "padj_hyper", "odds_ratio_log2_hyper", "upper_hyper", "lower_hyper")],
                           odds_ratio_results_hypo[,c("class_hypo", "padj_hypo", "odds_ratio_log2_hypo", "upper_hypo", "lower_hypo")])

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
