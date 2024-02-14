#'generate transposable element region metadata from regions of interest given provided annotation files
#'
#' @param regions regions of interest
#' @param conserved_annotation annotation file
#' @return Function returns region metadata of interest
#' @export human_conserved_region_metaData_generation

human_conserved_region_metaData_generation <- function(regions,
                                       conserved_annotation
) {
  chromosome <- unique(regions$chr)
  if (length(chromosome)>1){stop("Run this function in parallel or one chromosome at a time.")}

  i = 1
  for (region in rownames(regions)) {

    ph <- stringr::str_split(region, stringr::fixed("_"))
    start <-  as.numeric(ph[[1]][2])
    end <-  as.numeric(ph[[1]][3])

    # this is
    rangesA <- IRanges::IRanges(start, end)
    rangesB <- IRanges::IRanges(conserved_annotation$start, conserved_annotation$end)

    ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
    hit <- conserved_annotation[ov,]

    if (nrow(hit) > 0){conserved_region = 1} else {conserved_region = 0}

    if (i == 1){
      region_metaData <- data.frame("region" = region,"conserved_region"=conserved_region)
      i = 2
    } else {
      region_metaData <- rbind(region_metaData, data.frame("region" = region,"conserved_region"=conserved_region))
    }
  }
  return(region_metaData)
}
