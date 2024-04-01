#'generate transposable element region metadata from regions of interest given provided annotation files
#'
#' @param regions regions of interest
#' @param transposons_annotation TE annotation file
#' @param sw_score_cutoff the minimum sw score that will be queried here
#' @param expanded_search_if_nec expanding the search if necessary for region overlaps
#' @param return_classes if TRUE, returns if the region overlaps with transposon classes (LINE1, LINE2, SINE_tRNA, etc.)
#' @return Function returns region metadata of interest
#' @export TE_region_metaData_generation

TE_region_metaData_generation <- function(regions,
                                       transposons_annotation,
                                       sw_score_cutoff = 225,
                                       expanded_search_if_nec = 50,
                                       return_classes = FALSE
) {
  chromosome <- unique(regions$chr)
  if (length(chromosome)>1){stop("Run this function in parallel or one chromosome at a time.")}

  colnames(transposons_annotation) <- c('bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns'	,
                           'seqnames', 'start', 'end', 'genoLeft', 'strand'	,
                           'repName', 'class', 'repFamily', 'repStart', 'repEnd',	'repLeft',	'id')
  transposons_annotation <- transposons_annotation[transposons_annotation$seqnames == chromosome,]
  transposons_annotation <- transposons_annotation[transposons_annotation$swScore > sw_score_cutoff,]

  transposons_annotation$length <- transposons_annotation$end - transposons_annotation$start

  transposons_annotation$LINE_id <- 0
  transposons_annotation$SINE_id <- 0
  transposons_annotation$LINE_id[transposons_annotation$class == "LINE"] <- paste0(transposons_annotation$class[transposons_annotation$class == "LINE"], "_",
                                           transposons_annotation$repFamily[transposons_annotation$class == "LINE"], "_", transposons_annotation$repName[transposons_annotation$class == "LINE"])
  transposons_annotation$SINE_id[transposons_annotation$class == "SINE"] <- paste0(transposons_annotation$class[transposons_annotation$class == "SINE"], "_",
                                           transposons_annotation$repFamily[transposons_annotation$class == "SINE"], "_", transposons_annotation$repName[transposons_annotation$class == "SINE"])
  transposons_annotation$DNA_id[transposons_annotation$class == "DNA"] <- paste0(transposons_annotation$class[transposons_annotation$class == "DNA"], "_",
                                                                                   transposons_annotation$repFamily[transposons_annotation$class == "DNA"], "_", transposons_annotation$repName[transposons_annotation$class == "DNA"])
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
      rangesA <- IRanges::IRanges(start-expanded_search_if_nec, end+expanded_search_if_nec)
      ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
      hit <- transposons_annotation[ov,]
    }

    length_ = 0

    if ("DNA" %in% hit$class){dna = 1; length_ <- paste(unique(hit$length[hit$DNA_id!=0]), collapse = " & "); dna_id <- unique(hit$DNA_id[hit$DNA_id !=0])} else {dna = 0; dna_id = 0}
    if ("LTE" %in% hit$class){ltr = 1} else {ltr = 0}
    if ("rRNA" %in% hit$class){rRNA = 1} else {rRNA = 0}
    if ("tRNA" %in% hit$class){tRNA = 1} else {tRNA = 0}
    if ("Satellite" %in% hit$class){Satellite = 1} else {Satellite = 0}
    if ("Simple_repeat" %in% hit$class){Simple_repeat = 1} else {Simple_repeat = 0}
    if ("Low_complexity" %in% hit$class){Low_complexity = 1} else {Low_complexity = 0}
    if ("Unknown" %in% hit$class){Unknown = 1} else {Unknown = 0}
    if ("LINE" %in% hit$class){line = 1; length_ <- paste(unique(hit$length[hit$LINE_id!=0]), collapse = " & "); line_id <- unique(hit$LINE_id[hit$LINE_id!=0])} else {line = 0; line_id <- 0}
    if ("SINE" %in% hit$class){sine = 1; length_ <- paste(unique(hit$length[hit$SINE_id!=0]), collapse = " & "); sine_id <- unique(hit$SINE_id[hit$SINE_id!=0])} else {sine = 0; sine_id <- 0}

    if (grepl("L1_Ca", paste(line_id, collapse = " & ")) | grepl("L1M", paste(line_id, collapse = " & ")) ){l1_young <- 1}else{l1_young<-0}
    if (grepl("L1", paste(line_id, collapse = " & "))){l1 <- 1}else{l1<-0}
    if (grepl("L2", paste(line_id, collapse = " & "))){l2 <- 1}else{l2<-0}

    if (i == 1){
      region_metaData <- data.frame("region" = region,
                                    "DNA_transposon" = dna,
                                    "LTR" = ltr,
                                    "LINE" = line,
                                    "LINE1young" = l1_young,
                                    "LINE1" = l1,
                                    "LINE2" = l2,
                                    "SINE" = sine,
                                    "LINE_id" = paste(line_id, collapse = " & "),
                                    "SINE_id" = paste(sine_id, collapse = " & "),
                                    "DNA_id" = paste(dna_id, collapse = " & "),
                                    "length_LINESINEDNA" = length_,
                                    "rRNA" = paste(rRNA, collapse = " & "),
                                    "tRNA" = paste(tRNA, collapse = " & "),
                                    "Satellite" = paste(Satellite, collapse = " & "),
                                    "Simple_repeat" = paste(Simple_repeat, collapse = " & "),
                                    "Low_complexity" = paste(Low_complexity, collapse = " & "),
                                    "Unknown" = paste(Unknown, collapse = " & ")
                                    )

      i = 2
    } else {
      region_metaData <- rbind(region_metaData, data.frame("region" = region,
                                                           "DNA_transposon"= dna,
                                                           "LTR" = ltr,
                                                           "LINE" = line,
                                                           "LINE1young" = l1_young,
                                                           "LINE1" = l1,
                                                           "LINE2" = l2,
                                                           "SINE" = sine,
                                                           "LINE_id" = paste(line_id, collapse = " & "),
                                                           "SINE_id" = paste(sine_id, collapse = " & "),
                                                           "DNA_id" = paste(dna_id, collapse = " & "),
                                                           "length_DNALINESINE" = length_,
                                                           "rRNA" = paste(rRNA, collapse = " & "),
                                                           "tRNA" = paste(tRNA, collapse = " & "),
                                                           "Satellite" = paste(Satellite, collapse = " & "),
                                                           "Simple_repeat" = paste(Simple_repeat, collapse = " & "),
                                                           "Low_complexity" = paste(Low_complexity, collapse = " & "),
                                                           "Unknown" = paste(Unknown, collapse = " & ")
                                                           )
                               )
    }
  }

  region_metaData$Retrotransposon <- 0
  region_metaData$Retrotransposon[region_metaData$LINE == 1 | region_metaData$SINE == 1 | region_metaData$LTR == 1] <- 1

  region_metaData$TE <- 0
  region_metaData$TE[region_metaData$Retrotransposon == 1 | region_metaData$DNA_transposon == 1] <- 1

  return(region_metaData)
}
