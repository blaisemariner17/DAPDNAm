#'generate region metadata from regions of interest given provided annotation files
#'
#' @param regions regions of interest
#' @param transposons_annotation TE annotation file
#' @param sw_score_cutoff the minimum sw score that will be queried here
#' @return Function returns region metadata of interest
#' @export TE_region_metaData_generation

TE_region_metaData_generation <- function(regions,
                                       transposons_annotation,
                                       sw_score_cutoff = 0
) {
  chromosome <- unique(regions$chr)
  if (length(chromosome)>1){stop("Run this function in parallel or one chromosome at a time.")}

  colnames(transposons_annotation) <- c('bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns'	,
                           'seqnames', 'start', 'end', 'genoLeft', 'strand'	,
                           'repName', 'class', 'repFamily', 'repStart', 'repEnd',	'repLeft',	'id')
  transposons_annotation <- transposons_annotation[transposons_annotation$seqnames == chromosome,]
  transposons_annotation <- transposons_annotation[transposons_annotation$swScore > sw_score_cutoff,]
  transposons_annotation <- transposons_annotation[transposons_annotation$class %in% c("LINE", "SINE", "LTR",
                                                "Satellite", "tRNA",
                                                "snRNA", "rRNA", "scRNA", "srpRNA",
                                                #added 23-12-29
                                                "Simple_repeat", "Low_complexity", "DNA", "RC"),]
  transposons_annotation$LINE_id <- 0
  transposons_annotation$SINE_id <- 0
  transposons_annotation$LINE_id[transposons_annotation$class == "LINE"] <- paste0(transposons_annotation$class[transposons_annotation$class == "LINE"], "_",
                                           transposons_annotation$repFamily[transposons_annotation$class == "LINE"], "_",
                                           transposons_annotation$id[transposons_annotation$class == "LINE"])
  transposons_annotation$SINE_id[transposons_annotation$class == "SINE"] <- paste0(transposons_annotation$class[transposons_annotation$class == "SINE"], "_",
                                           transposons_annotation$repFamily[transposons_annotation$class == "SINE"], "_",
                                           transposons_annotation$id[transposons_annotation$class == "SINE"])
  i = 1
  for (region in rownames(regions)) {

    ph <- stringr::str_split(region, stringr::fixed("_"))
    start <-  as.numeric(ph[[1]][2])
    end <-  as.numeric(ph[[1]][3])

    # this is
    rangesA <- IRanges::IRanges(start, end)
    rangesB <- IRanges::IRanges(transposons_annotation$start, transposons_annotation$end)

    ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
    hit <- transposons_annotation[ov,]
    #expand the search of the region to assess if there is a nearby annotation in an effort to improve our region annotation
    if (nrow(hit) == 0){
      rangesA <- IRanges::IRanges(start-1000, end+1000)
      ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
      hit <- transposons_annotation[ov,]
    }
    sine_id <- 0
    line_id <- 0

    if ("LINE" %in% hit$class){line = 1; line_id <- unique(hit$LINE_id[hit$LINE_id!=0])} else {line = 0}
    if ("SINE" %in% hit$class){sine = 1; sine_id <- unique(hit$SINE_id[hit$SINE_id!=0])} else {sine = 0}

    if (i == 1){
      region_metaData <- data.frame("region" = region,"LINE" = line,"SINE" = sine,"LINE_id" = paste(line_id, collapse = " & "),"SINE_id" = paste(sine_id, collapse = " & "))
      i = 2
    } else {
      region_metaData <- rbind(region_metaData, data.frame("region" = region,"LINE" = line,"SINE" = sine,"LINE_id" = paste(line_id, collapse = " & "),"SINE_id" = paste(sine_id, collapse = " & ")))
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
