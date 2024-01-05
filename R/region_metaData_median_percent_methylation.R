#'returns a list of plots
#'
#' @param coverage_all_chr region coverage
#' @param methylation_all_chr region methylation
#' @param region_metaData region metaData
#' @param col_oi columns of interest of the region metaData
#' @param theme_blaise ggplot theme for plotting
#' @return Function returns a list of plots
#' @export region_metaData_generation

region_metaData_median_percent_methylation <- function(coverage_all_chr,
                                                       methylation_all_chr,
                                                       region_metaData,
                                                       col_oi = c( "Promoter","exon","intron","upstream_utr","downstram_utr",
                                                                    "CpG_shore","CpG_island",#"CpG_shelf",
                                                                    "LINE","SINE",
                                                                    #promoters
                                                                    "ActiveTSS","WeakTSS","FlankingActiveTSS1","FlankingActiveTSS2",
                                                                    #enhancers
                                                                    "ActiveStrongEnhancer","ActiveWeakEnhancer", "ActivePoisedEnhancer",
                                                                    ##polycomb-repression
                                                                    "BivalentTssEnh","RepressedPolycomb",
                                                                    #Heterochromatin
                                                                    "Repressed","ZNFgenesRepeats","QuiescenatLow","Heterochromatin"
                                                       ),
                                                       theme_blaise = theme(axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 24, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      axis.text = element_text(family = "sans", size = 18, color="black"),       axis.title = element_text(family = "sans", size = 20, color="black"),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

) {

  percent_methylation <- methylation_all_chr / coverage_all_chr

  plots <- list()
  for (col in col_oi){#colnames(region_metaData)){
    if (col == "gene_bool"){
      next
    }
    if (T){#1 %in% region_metaData[,c(paste(col))]){
      print(col)
      region_metaData_oi <- region_metaData[region_metaData[,c(paste(col))] == 1,]
      percent_methylation_oi <- percent_methylation[rownames(percent_methylation) %in% region_metaData_oi$region,]
      for_ggplot <- data.frame(rowMedians(percent_methylation_oi, na.rm = T))
      colnames(for_ggplot) <- c("data")
      plot_ <- ggplot(for_ggplot, aes(x = data)) +
        geom_density(fill = "lightgrey") +
        xlab("")+
        ylab(paste0(col, " (n = ", nrow(for_ggplot), ")"))+
        ggtitle("") +
        theme_blaise +
        xlim(c(0,1))+
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.title.y = element_text(angle = 0, hjust = -1, vjust = 0.5),
              axis.title.x = element_blank(), plot.title = element_blank(),
              plot.margin = margin(0, 0, 0, 0, "cm")
        )
      plot_
      if (col == "QuiescenatLow") {
        plot_ <- plot_ + xlab("Median fraction methylated") + theme(axis.text.x=element_text(family = "sans", size = 18, color="black"),
                                                                    axis.line.x = element_line(colour = "black", linewidth = 1),
                                                                    axis.title.x = element_text(family = "sans", size = 24, hjust = 0.5, color="black")
        )
      } else {
        plot_ <- plot_ + theme(axis.line.x = element_blank())
      }

      plots[[col]] <- plot_
    }
  }

  return(plots)
}
