#'generate region metadata from regions of interest given provided annotation files
#'
#' @param regions regions of interest
#' @param gtf gene annotation file
#' @param gtf_trans transposon annotation file
#' @param gtf_cpgisl CpG island annotation file
#' @param chromatin_state_bed_file bed file of the chromatin state annotation
#' @param n.cores numer of cores
#' @return Function returns region metadata of interest
#' @export region_metaData_generation_old

region_metaData_generation_old <- function(regions,
                                       n.cores,
                                       gtf_trans = gtf_trans,
                                       gtf = rtracklayer::import('../../GENOME-ANNOTATION-FILE/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf'),
                                       gtf_cpgisl = gtf_cpgisl,
                                       chromatin_state_bed_file = bed_file
) {

  #### First, let's generate cpg shalf and shore annotations
  gtf_cpgisl$class <- "CpG_island"
  colnames(gtf_cpgisl) <- c("idk", "seqnames", "start", "end", "id", "length", "Num_CpG",
                            "idk","idk","idk","idk", "class")
  #add the cpg shores, which are the regions flanking the islands by as much as 2kb
  gtf_cpgisl$shoreup_start <- gtf_cpgisl$start - 2001
  gtf_cpgisl$shoreup_end <- gtf_cpgisl$start-1
  gtf_cpgisl$shelfup_start <- gtf_cpgisl$shoreup_start - 2001
  gtf_cpgisl$shelfup_end <- gtf_cpgisl$shoreup_start - 1

  gtf_cpgisl$shoredown_start <- gtf_cpgisl$end+1
  gtf_cpgisl$shoredown_end <- gtf_cpgisl$end + 2001
  gtf_cpgisl$shelfdown_start <- gtf_cpgisl$shoredown_end + 1
  gtf_cpgisl$shelfdown_end <- gtf_cpgisl$shoredown_end + 2001

  gtf_cpgshoreup <- gtf_cpgisl[,c("seqnames", "shoreup_start", "shoreup_end", "id")]
  gtf_cpgshoredown <- gtf_cpgisl[,c("seqnames", "shoredown_start", "shoredown_end", "id")]
  gtf_cpgshelfup <- gtf_cpgisl[,c("seqnames", "shelfup_start", "shelfup_end", "id")]
  gtf_cpgshelfdown <- gtf_cpgisl[,c("seqnames", "shelfdown_start", "shelfdown_end", "id")]

  colnames(gtf_cpgshoreup) <- c("seqnames", "start", "end","id")
  colnames(gtf_cpgshoredown) <- c("seqnames", "start", "end","id")
  colnames(gtf_cpgshelfup) <- c("seqnames", "start", "end","id")
  colnames(gtf_cpgshelfdown) <- c("seqnames", "start", "end","id")

  gtf_cpgshores <- rbind(gtf_cpgshoreup, gtf_cpgshoredown)
  gtf_cpgshelves <- rbind(gtf_cpgshelfup, gtf_cpgshelfdown)

  #determine if there is a cpg island in a shore and remove the shore if so
  counter = 1
  for (chr in c(1:38,"X")){
    chromosome = paste0("chr",chr)
    cpg_island_chroi <- as.data.frame(gtf_cpgisl[gtf_cpgisl$seqnames == chromosome,])
    cpg_shores_chroi <- as.data.frame(gtf_cpgshores[gtf_cpgshores$seqnames == chromosome,])

    rangesA <- split(IRanges(cpg_island_chroi$start, cpg_island_chroi$end), chromosome)
    rangesB <- split(IRanges(cpg_shores_chroi$start, cpg_shores_chroi$end), chromosome)

    #which regionsB overlap w no regionA regions
    not_ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type = 'any')>0

    cpg_shores_chroi <- cpg_shores_chroi[not_ov[[1]],]
    cpg_shores_chroi$class <- "CpG_shore"

    if (counter == 1) {
      res <- cpg_shores_chroi
      counter = 2
    } else {
      res <- rbind(res, cpg_shores_chroi)
    }
  }
  gtf_cpgshores <- res

  #determine if there is a cpg shore in a shelf and remove the shelf if so
  counter = 1
  for (chr in c(1:38,"X")){
    chromosome = paste0("chr",chr)
    cpg_shores_chroi <- as.data.frame(gtf_cpgshores[gtf_cpgshores$seqnames == chromosome,])
    cpg_shelf_chroi <- as.data.frame(gtf_cpgshelves[gtf_cpgshelves$seqnames == chromosome,])

    rangesA <- split(IRanges(cpg_shores_chroi$start, cpg_shores_chroi$end), chromosome)
    rangesB <- split(IRanges(cpg_shelf_chroi$start, cpg_shelf_chroi$end), chromosome)

    #which regionsB overlap w no regionA regions
    not_ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type = 'any')>0

    cpg_shelf_chroi <- cpg_shelf_chroi[not_ov[[1]],]
    cpg_shelf_chroi$class <- "CpG_shelf"

    if (counter == 1) {
      res <- cpg_shelf_chroi
      counter = 2
    } else {
      res <- rbind(res, cpg_shelf_chroi)
    }
  }
  gtf_cpgshelves <- res

  #combine all the cpg annotation
  gtf_cpgisl <- gtf_cpgisl[,c("seqnames", "start", "end","id", "class")]
  gtf_cpgisl <- rbind(gtf_cpgisl, gtf_cpgshores, gtf_cpgshelves)

  #### now let's look at the annotations from the transposon annotation provided
  colnames(gtf_trans) <- c('bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns'	,
                           'seqnames', 'start', 'end', 'genoLeft', 'strand'	,
                           'repName', 'class', 'repFamily', 'repStart', 'repEnd',	'repLeft',	'id')

  gtf_trans$id <- paste0(gtf_trans$class, "_", gtf_trans$repFamily, "_", gtf_trans$id )

  gtf_trans <- gtf_trans[gtf_trans$class %in% c("LINE", "SINE", "LTR",
                                                "Satellite", "tRNA",
                                                "snRNA", "rRNA", "scRNA", "srpRNA",
                                                #added 23-12-29
                                                "Simple_repeat", "Low_complexity", "DNA", "RC"),]

  #### now let's look at the chromatin states
  chromatin_states <- rtracklayer::import(chromatin_state_bed_file, format = "bed")

  ### now we must build the promoter region annotation
  promoter_gtf <- data.frame(gtf[gtf$type == "gene",])

  promoter_gtf_neg <- promoter_gtf[promoter_gtf$strand == "-",]
  promoter_gtf_neg$promoter_start <- promoter_gtf_neg$end
  promoter_gtf_neg$promoter_end <- promoter_gtf_neg$end + 2000

  promoter_gtf_pos <- promoter_gtf[promoter_gtf$strand == "+",]
  promoter_gtf_pos$promoter_start <- promoter_gtf_pos$start - 2000
  promoter_gtf_pos$promoter_end <- promoter_gtf_pos$start

  promoter_gtf <- rbind(promoter_gtf_pos, promoter_gtf_neg)
  promoter_gtf$type <- "Promoter"
  promoter_gtf$length <- promoter_gtf$promoter_end - promoter_gtf$promoter_start
  promoter_gtf <- promoter_gtf[,c("seqnames", "promoter_start", "promoter_end", "type", "gene_id")]

  ####I want is the regions df for each chromosome... and the mclapply function does not like when I write out as csv so I guess I have to do it the way I know how to

  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)

  library_list <-  c(
    "ggplot2",
    "svglite",
    "ggExtra",
    "ggtext",
    "GenomicRanges",
    "GenomicFeatures",
    "tidyverse",
    "bsseq",
    "comethyl",
    "PQLseq",
    "foreach",
    "parallel",
    "Biostrings",
    "patchwork"
  )

  time_start <- Sys.time()
  region_metaData <- foreach(region = rownames(regions),
                             .packages = library_list,
                             .combine = "rbind") %dopar% {
                               # print(region)
                               # build the annotation for all the gtf files you have for the chromosome that this region is on...
                               ph <- str_split(region, fixed("_"))
                               chromosome <- ph[[1]][1]
                               start <-  as.numeric(ph[[1]][2])
                               end <-  as.numeric(ph[[1]][3])
                               chr_oi_gtf <- as.data.frame(gtf[seqnames(gtf) == chromosome,])
                               chr_oi_gtf <- chr_oi_gtf[,c("seqnames", "start", "end", "type", "gene_id")]

                               promoter_gtf_oi <- promoter_gtf[promoter_gtf$seqnames == chromosome,]

                               colnames(promoter_gtf_oi) <- c("seqnames", "start", "end", "class", "id")
                               colnames(gtf_cpgisl) <- c("seqnames", "start", "end", "class", "id")
                               gtf_trans <- gtf_trans[,c("seqnames", "start", "end", "class", "id")]
                               colnames(chr_oi_gtf) <- c("seqnames", "start", "end", "class", "id")

                               chr_oi_all <- rbind(chr_oi_gtf,
                                                   promoter_gtf_oi,
                                                   gtf_cpgisl[gtf_cpgisl$seqnames == chromosome,],
                                                   gtf_trans[gtf_trans$seqnames == chromosome,]
                               )

                               rangesA <- split(IRanges(start, end), chromosome)
                               rangesB <- split(IRanges(chr_oi_all$start, chr_oi_all$end), chromosome)

                               #which rangesB have an overlap with at least one rangesA?
                               ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
                               hit <- chr_oi_all[ov[[1]],]

                               gene_id <- 0
                               intron_gene <- "no annotation"
                               gene_bool <- 0
                               closest_distance_to_promoter <- 0
                               Promoter <- 0
                               Promoter_id <- 0
                               downstram_utr <- 0
                               upstream_utr <- 0
                               exon <- 0
                               intron <- 0
                               LINE <- 0
                               LINE_id <- 0
                               SINE <- 0
                               SINE_id <- 0
                               Cpg_shelf <- 0
                               CpG_shore <- 0
                               CpG_island <- 0
                               Chromatin_states1 <- 0
                               Chromatin_states2 <- 0
                               Chromatin_states3 <- 0
                               Chromatin_states4 <- 0
                               Chromatin_states5 <- 0
                               Chromatin_states6 <- 0
                               Chromatin_states7 <- 0
                               Chromatin_states8 <- 0
                               Chromatin_states9 <- 0
                               Chromatin_states10 <- 0
                               Chromatin_states11 <- 0
                               Chromatin_states12 <- 0
                               Chromatin_states13 <- 0

                               gene_gtf_for_distance <- chr_oi_gtf[chr_oi_gtf$class == "gene",]
                               rangesB <- split(IRanges(gene_gtf_for_distance$start, gene_gtf_for_distance$end), chromosome)
                               query <- GRanges("A", IRanges(c(start, end), width=1), strand="+")
                               subject <- GRanges("A", IRanges(gene_gtf_for_distance$start))
                               d <- data.frame(distanceToNearest(query, subject, ignore.strand=TRUE))
                               distance_to_nearest_gene_start <- min(d$distance)

                               for (id in unique(hit$id)) {
                                 # print(id)
                                 if ("CpG_island" %in% hit$id[hit$id == id]) {CpG_island <- 1}
                                 if ("CpG_shore" %in% hit$id[hit$id == id]) {CpG_shore <- 1}
                                 if ("CpG_shelf" %in% hit$id[hit$id == id]) {CpG_shelf <- 1}
                                 if ("Simple_repeat" %in% hit$class[hit$id == id]) {Simple_repeat <- 1}
                                 if (TRUE %in% grepl("LINE", hit$class[hit$id == id])){
                                   LINE <- 1
                                   LINE_id <-  paste(unique(hit$id[hit$class == "LINE"]), collapse = " & ")
                                 }
                                 if (TRUE %in% grepl("SINE", hit$class[hit$id == id])) {
                                   SINE <- 1
                                   SINE_id <-  paste(unique(hit$id[hit$class == "SINE"]), collapse = " & ")
                                 }
                                 if ("Promoter" %in% hit$class[hit$id == id]) {
                                   Promoter <- 1
                                   Promoter_id <- paste(unique(id), collapse = ' & ')
                                 }

                                 #### chromatin states

                                 # Categorize cpg_sites based on chromatin state map ranges and names
                                 chromatin_states_oi <- as.data.frame(chromatin_states[seqnames(chromatin_states) == chromosome,])
                                 rangesB <- split(IRanges(chromatin_states_oi$start, chromatin_states_oi$end), chromosome)
                                 chromatin_overlaps <- GenomicRanges::countOverlaps(rangesB, rangesA, type = "any")>0
                                 hit_chromatin <- chromatin_states_oi[chromatin_overlaps[[1]],]
                                 # hit_chromatin
                                 if (nrow(hit_chromatin) > 0) {
                                   if (1 %in% hit_chromatin$name) {Chromatin_states1 <- 1}
                                   if (2 %in% hit_chromatin$name) {Chromatin_states2 <- 1}
                                   if (3 %in% hit_chromatin$name) {Chromatin_states3 <- 1}
                                   if (4 %in% hit_chromatin$name) {Chromatin_states4 <- 1}
                                   if (5 %in% hit_chromatin$name) {Chromatin_states5 <- 1}
                                   if (6 %in% hit_chromatin$name) {Chromatin_states6 <- 1}
                                   if (7 %in% hit_chromatin$name) {Chromatin_states7 <- 1}
                                   if (8 %in% hit_chromatin$name) {Chromatin_states8 <- 1}
                                   if (9 %in% hit_chromatin$name) {Chromatin_states9 <- 1}
                                   if (10 %in% hit_chromatin$name) {Chromatin_states10 <- 1}
                                   if (11 %in% hit_chromatin$name) {Chromatin_states11 <- 1}
                                   if (12 %in% hit_chromatin$name) {Chromatin_states12 <- 1}
                                   if (13 %in% hit_chromatin$name) {Chromatin_states13 <- 1}
                                 }

                                 hit_intron_check <- hit[hit$id == id,]
                                 if ("gene" %in% hit_intron_check$class){
                                   gene_id <- paste(unique(hit_intron_check$id[hit_intron_check$class == "gene"]), collapse = " & ")
                                   gene_bool <- 1
                                   if ("exon" %in% chr_oi_all$class[chr_oi_all$id == id]) {
                                     if (sum(grepl("exon", chr_oi_all$class[chr_oi_all$id == id])) == 1){
                                       intron_gene <- "no"
                                     } else {intron_gene <- "yes"}
                                   } else {intron_gene <- "no"}
                                   if ("exon" %in% hit_intron_check$class){
                                     exon <- 1
                                     next
                                   }
                                   if ("start_codon" %in% chr_oi_all$class[chr_oi_all$id == id]) {
                                     start_codon <- chr_oi_all$start[chr_oi_all$id == id & chr_oi_all$class == "start_codon"][1]
                                     stop_codon <- chr_oi_all$start[chr_oi_all$id == id & chr_oi_all$class == "stop_codon"][1]
                                     # what region of the gene is it in if it is not in the exon? have to account for the instances the gene is on the negative or the positive strand here and if the annotation file even documents exon or intron
                                     if (start_codon < stop_codon) {
                                       if (end < start_codon) {upstream_utr <- 1} else if (start > stop_codon) {downstram_utr <- 1} else if ("exon" %in% chr_oi_all$class[chr_oi_all$id == id]) {intron <- 1}
                                     } else if (start_codon > stop_codon) {
                                       if (start > start_codon) {upstream_utr <- 1} else if (end < stop_codon) {downstram_utr <- 1} else if ("exon" %in% chr_oi_all$class[chr_oi_all$id == id]) {intron <- 1}
                                     }
                                   }
                                 }
                               }
                               metaData_ <- data.frame("region" = region,
                                                       "gene_id" = gene_id,
                                                       "intron_gene" = intron_gene,
                                                       "distance_nearest_gene_start" = distance_to_nearest_gene_start,
                                                       "closest_distance_to_promoter" = closest_distance_to_promoter,
                                                       "gene_bool" = gene_bool,
                                                       "Promoter" = Promoter,
                                                       "Promoter_id" = Promoter_id,
                                                       "exon" = exon,
                                                       "intron" = intron,
                                                       "upstream_utr" = upstream_utr,
                                                       "downstram_utr" = downstram_utr,
                                                       "LINE" = LINE,
                                                       "LINE_id" = LINE_id,
                                                       "SINE" = SINE,
                                                       "SINE_id" = SINE_id,
                                                       "CpG_shelf" = Cpg_shelf,
                                                       "CpG_shore" = CpG_shore,
                                                       "CpG_island" = CpG_island,
                                                       "ActiveTSS" = Chromatin_states1,
                                                       "WeakTSS" = Chromatin_states2,
                                                       "FlankingActiveTSS1" = Chromatin_states3,
                                                       "FlankingActiveTSS2" = Chromatin_states4,
                                                       "ActiveStrongEnhancer" = Chromatin_states5,
                                                       "ActiveWeakEnhancer" = Chromatin_states6,
                                                       "ActivePoisedEnhancer" = Chromatin_states7,
                                                       "BivalentTssEnh" = Chromatin_states8,
                                                       "RepressedPolycomb" = Chromatin_states9,
                                                       "Repressed" = Chromatin_states10,
                                                       "ZNFgenesRepeats" = Chromatin_states11,
                                                       "Heterochromatin" = Chromatin_states12,
                                                       "QuiescenatLow" = Chromatin_states13

                               )
                               return(metaData_)
                             }
  print(Sys.time() - time_start)
  parallel::stopCluster(cl = my.cluster)

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

  return(region_metaData)
}
