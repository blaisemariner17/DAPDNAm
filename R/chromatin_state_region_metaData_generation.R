#'generate chromatin region metadata from regions of interest given provided annotation file
#'
#' @param regions regions of interest
#' @param chromatin_states_annotation chromatin states annotation
#' @param n.cores number of cores
#' @return Function returns region metadata of interest
#' @export chromatin_state_region_metaData_generation

chromatin_state_region_metaData_generation <- function(regions,
                                                       chromatin_states_annotation
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

  return(region_metaData)
}