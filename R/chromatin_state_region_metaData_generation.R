#'generate region metadata from regions of interest given provided annotation files
#'
#' @param regions regions of interest
#' @param chromatin_states_annotation chromatin states annotation
#' @param n.cores numer of cores
#' @return Function returns region metadata of interest
#' @export chromatin_state_region_metaData_generation

chromatin_state_region_metaData_generation <- function(regions,
                                       genome_gene_annotation,
                                       cpgisland_annotation
) {
  chromosome <- unique(regions$chr)
  if (length(chromosome)>1){stop("Run this function in parallel or one chromosome at a time.")}
  
  #chromatin state annotation
  chromatin_states_annotation <- chromatin_states_annotation[chromatin_states_annotation$seqnames == chromosome,]
  
  i = 1
  for (region in rownames(regions)) {
    
    ph <- stringr::str_split(region, stringr::fixed("_"))
    start <-  as.numeric(ph[[1]][2])
    end <-  as.numeric(ph[[1]][3])
    
    rangesA <- IRanges::IRanges(start, end)
    rangesB <- IRanges::IRanges(chromatin_states_annotation$start, chromatin_states_annotation$end)
    
    ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
    hit <- chromatin_states_annotation[ov,]
    
    if ("Chromatin_state_promoter" %in% hit$class){Chromatin_state_promoter = 1} else {Chromatin_state_promoter = 0}
    if ("Chromatin_state_enhancer" %in% hit$class){Chromatin_state_enhancer = 1} else {Chromatin_state_enhancer = 0}
    if ("Chromatin_state_polycomb" %in% hit$class){Chromatin_state_polycomb = 1} else {Chromatin_state_polycomb = 0}
    if ("Chromatin_state_heterochromatin" %in% hit$class){Chromatin_state_heterochromatin = 1} else {Chromatin_state_heterochromatin = 0}
    if ("Chromatin_state_quies" %in% hit$class){Chromatin_state_quies = 1} else {Chromatin_state_quies = 0}
    
    if (i == 1){
      region_metaData <- data.frame("region" = region, "Chromatin_state_promoter" = Chromatin_state_promoter,
                                    "Chromatin_state_enhancer" = Chromatin_state_enhancer, "Chromatin_state_polycomb" = Chromatin_state_polycomb,
                                    "Chromatin_state_heterochromatin" = Chromatin_state_heterochromatin, "Chromatin_state_quies" = Chromatin_state_quies
                                    )
      i = 2
    } else {
      region_metaData <- rbind(region_metaData, data.frame("region" = region, "Chromatin_state_promoter" = Chromatin_state_promoter,
                                                           "Chromatin_state_enhancer" = Chromatin_state_enhancer, "Chromatin_state_polycomb" = Chromatin_state_polycomb,
                                                           "Chromatin_state_heterochromatin" = Chromatin_state_heterochromatin, "Chromatin_state_quies" = Chromatin_state_quies
      ))
    }
  }
  
  region_metaData$LINE1 <- 0
  region_metaData$LINE1[grepl("L1", region_metaData$LINE_id)] <- 1
  region_metaData$LINE2 <- 0
  region_metaData$LINE2[grepl("L2", region_metaData$LINE_id)] <- 1
  region_metaData$LINE3 <- 0
  region_metaData$LINE3[grepl("L3", region_metaData$LINE_id)] <- 1
  region_metaData$LINE_CR1 <- 0
  region_metaData$LINE_CR1[grepl("CR1", region_metaData$LINE_id)] <- 1
  
  region_metaData$SINEC <- 0
  region_metaData$SINEC[grepl("SINEC", region_metaData$SINE_id)] <- 1
  region_metaData$SINE_MIR <- 0
  region_metaData$SINE_MIR[grepl("MIR", region_metaData$SINE_id)] <- 1
  region_metaData$SINE_tRNA <- 0
  region_metaData$SINE_tRNA[grepl("tRNA", region_metaData$SINE_id)] <- 1
  
  region_metaData$TE <- 0
  region_metaData$TE[region_metaData$LINE == 1 | region_metaData$SINE == 1] <- 1
  
  return(region_metaData)
  }
}