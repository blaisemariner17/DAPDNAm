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
  chromatin_states_annotation <- chromatin_states_annotation[chromatin_states_annotation@seqnames == chromosome,]

  i = 1
  for (region in rownames(regions)) {

    ph <- stringr::str_split(region, stringr::fixed("_"))
    start <-  as.numeric(ph[[1]][2])
    end <-  as.numeric(ph[[1]][3])

    rangesA <- IRanges::IRanges(start, end)
    rangesB <- IRanges::IRanges(chromatin_states_annotation@ranges)

    ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
    hit <- chromatin_states_annotation[ov,]

    if ("ChrSt_promoter" %in% hit$name){ChrSt_promoter = 1} else {ChrSt_promoter = 0}
    if ("ChrSt_enhancer" %in% hit$name){ChrSt_enhancer = 1} else {ChrSt_enhancer = 0}
    if ("ChrSt_polycomb" %in% hit$name){ChrSt_polycomb = 1} else {ChrSt_polycomb = 0}
    if ("ChrSt_heterochromatin" %in% hit$name){ChrSt_heterochromatin = 1} else {ChrSt_heterochromatin = 0}
    if ("ChrSt_quies" %in% hit$name){ChrSt_quies = 1} else {ChrSt_quies = 0}

    if (i == 1){
      region_metaData <- data.frame("region" = region, "ChrSt_promoter" = ChrSt_promoter,
                                    "ChrSt_enhancer" = ChrSt_enhancer, "ChrSt_polycomb" = ChrSt_polycomb,
                                    "ChrSt_heterochromatin" = ChrSt_heterochromatin, "ChrSt_quies" = ChrSt_quies
                                    )
      i = 2
    } else {
      region_metaData <- rbind(region_metaData, data.frame("region" = region, "ChrSt_promoter" = ChrSt_promoter,
                                                           "ChrSt_enhancer" = ChrSt_enhancer, "ChrSt_polycomb" = ChrSt_polycomb,
                                                           "ChrSt_heterochromatin" = ChrSt_heterochromatin, "ChrSt_quies" = ChrSt_quies
      ))
    }
  }

  return(region_metaData)
}
