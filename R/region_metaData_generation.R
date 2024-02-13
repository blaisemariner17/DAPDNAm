#'generate region metadata from regions of interest given provided annotation files
#'
#' @param regions regions of interest
#' @param genome_gene_annotation gene annotation file
#' @param cpgisland_annotation CpG island annotation file
#' @param n.cores numer of cores
#' @return Function returns region metadata of interest
#' @export region_metaData_generation

region_metaData_generation <- function(regions,
                                       genome_gene_annotation,
                                       cpgisland_annotation
) {

  chromosome <- unique(regions$chr)
  if (length(chromosome)>1){stop("Run this function in parallel or one chromosome at a time.")}

  #building the promoter annotation
  promoter_annotation <- data.frame(genome_gene_annotation[genome_gene_annotation$type == "gene",])
  promoter_annotation <- data.frame(promoter_annotation[promoter_annotation$seqnames == chromosome,])

  promoter_annotation_neg <- promoter_annotation[promoter_annotation$strand == "-",]
  promoter_annotation_neg$promoter_start <- promoter_annotation_neg$end
  promoter_annotation_neg$promoter_end <- promoter_annotation_neg$end + 2000

  promoter_annotation_pos <- promoter_annotation[promoter_annotation$strand == "+",]
  promoter_annotation_pos$promoter_start <- promoter_annotation_pos$start - 2000
  promoter_annotation_pos$promoter_end <- promoter_annotation_pos$start

  promoter_annotation <- rbind(promoter_annotation_pos, promoter_annotation_neg)
  promoter_annotation$type <- "Promoter"
  promoter_annotation$length <- promoter_annotation$promoter_end - promoter_annotation$promoter_start
  promoter_annotation <- promoter_annotation[,c("seqnames", "promoter_start", "promoter_end", "type", "gene_id")]
  colnames(promoter_annotation) <- c("seqnames", "start", "end", "class", "id")

  #### Now, the cpg shelf and shore annotations
  cpgisland_annotation$class <- "CpG_island"
  colnames(cpgisland_annotation) <- c("idk", "seqnames", "start", "end", "id", "length", "Num_CpG",
                            "idk","idk","idk","idk", "class")
  cpgisland_annotation <- cpgisland_annotation[cpgisland_annotation$seqnames == chromosome,]
  #add the cpg shores, which are the regions flanking the islands by as much as 2kb
  cpgisland_annotation$shoreup_start <- cpgisland_annotation$start - 2001
  cpgisland_annotation$shoreup_end <- cpgisland_annotation$start-1
  cpgisland_annotation$shelfup_start <- cpgisland_annotation$shoreup_start - 2001
  cpgisland_annotation$shelfup_end <- cpgisland_annotation$shoreup_start - 1

  cpgisland_annotation$shoredown_start <- cpgisland_annotation$end+1
  cpgisland_annotation$shoredown_end <- cpgisland_annotation$end + 2001
  cpgisland_annotation$shelfdown_start <- cpgisland_annotation$shoredown_end + 1
  cpgisland_annotation$shelfdown_end <- cpgisland_annotation$shoredown_end + 2001

  gtf_cpgshoreup <- cpgisland_annotation[,c("seqnames", "shoreup_start", "shoreup_end", "id")]
  gtf_cpgshoredown <- cpgisland_annotation[,c("seqnames", "shoredown_start", "shoredown_end", "id")]
  gtf_cpgshelfup <- cpgisland_annotation[,c("seqnames", "shelfup_start", "shelfup_end", "id")]
  gtf_cpgshelfdown <- cpgisland_annotation[,c("seqnames", "shelfdown_start", "shelfdown_end", "id")]

  colnames(gtf_cpgshoreup) <- c("seqnames", "start", "end","id")
  colnames(gtf_cpgshoredown) <- c("seqnames", "start", "end","id")
  colnames(gtf_cpgshelfup) <- c("seqnames", "start", "end","id")
  colnames(gtf_cpgshelfdown) <- c("seqnames", "start", "end","id")

  gtf_cpgshores <- rbind(gtf_cpgshoreup, gtf_cpgshoredown)
  gtf_cpgshelves <- rbind(gtf_cpgshelfup, gtf_cpgshelfdown)

  #determine if there is a cpg island in a shore and remove the shore if so
  rangesA <- IRanges::IRanges(cpgisland_annotation$start, cpgisland_annotation$end)
  rangesB <- IRanges::IRanges(gtf_cpgshores$start, gtf_cpgshores$end)

  #which regionsB overlap w no regionA regions
  not_ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type = 'any')==0

  gtf_cpgshores <- gtf_cpgshores[not_ov,]
  gtf_cpgshores$class <- "CpG_shore"

  #determine if there is a cpg shore in a shelf and remove the shelf if so

  rangesA <- IRanges::IRanges(gtf_cpgshores$start, gtf_cpgshores$end)
  rangesB <- IRanges::IRanges(gtf_cpgshelves$start, gtf_cpgshelves$end)

  #which regionsB overlap w no regionA regions
  not_ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type = 'any')==0

  gtf_cpgshelves <- gtf_cpgshelves[not_ov,]
  gtf_cpgshelves$class <- "CpG_shelf"

  #combine all the cpg annotations
  cpgisland_annotation <- cpgisland_annotation[,c("seqnames", "start", "end","id", "class")]
  cpgisland_annotation <- rbind(cpgisland_annotation, gtf_cpgshores, gtf_cpgshelves)

  #gene body annotation
  gene_body_annotation <- as.data.frame(genome_gene_annotation[GenomeInfoDb::seqnames(genome_gene_annotation) == chromosome,])
  gene_body_annotation <- gene_body_annotation[,c("seqnames", "start", "end", "type", "gene_id")]
  colnames(gene_body_annotation) <- c("seqnames", "start", "end", "class", "id")

  all_compiled_annotations <- rbind(gene_body_annotation,
                      promoter_annotation,
                      cpgisland_annotation
  )

  i = 1
  for (region in rownames(regions)){
    ph <- stringr::str_split(region, stringr::fixed("_"))
    start <-  as.numeric(ph[[1]][2])
    end <-  as.numeric(ph[[1]][3])

    rangesA <- IRanges::IRanges(start, end)
    rangesB <- IRanges::IRanges(all_compiled_annotations$start, all_compiled_annotations$end)

    #which regionsB overlap w regionA
    ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
    hit <- all_compiled_annotations[ov,]

    #expand our search, if necessary
    if (nrow(hit) == 0){
      rangesA <- IRanges::IRanges(start-1000, end+1000)
      ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
      hit <- all_compiled_annotations[ov,]
    }
    if (nrow(hit) == 0){
      rangesA <- IRanges::IRanges(start-2000, end+2000)
      ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type="any")>0
      hit <- all_compiled_annotations[ov,]
    }

    gene_id <- 0
    gene_bool <- 0
    exon <- 0
    intron <- 0
    intron_gene <- 0
    upstream_utr <- 0
    downstream_utr <- 0
    Promoter <- 0
    Promoter_id <- 0

    if ("CpG_island" %in% hit$class) {CpG_island <- 1} else {CpG_island <- 0}
    if ("CpG_shore" %in% hit$class) {CpG_shore <- 1} else {CpG_shore <- 0}
    if ("CpG_shelf" %in% hit$class) {CpG_shelf <- 1} else {CpG_shelf <- 0}
    if ("Simple_repeat" %in% hit$class) {Simple_repeat <- 1} else {Simple_repeat <- 0}
    if ("Promoter" %in% hit$class) {
      Promoter <- 1
      Promoter_id <- paste(unique(hit$id[hit$class == "Promoter"]), collapse = ' & ')
    }
    if ("gene" %in% hit$class) {
      gene_bool <- 1
      gene_id <- paste(unique(hit$id[hit$class == "gene"]), collapse = ' & ')
    }

    if(nrow(hit[hit$class == 'gene',]) > 0){
      for (id in unique(hit$id)[hit$class == 'gene']){
        hit_intron_check <- hit[hit$id == id,]
        if ("exon" %in% all_compiled_annotations$class[all_compiled_annotations$id == id]) {
          if (sum(grepl("exon", all_compiled_annotations$class[all_compiled_annotations$id == id])) == 1){
            intron_gene <- "no"
          } else {intron_gene <- "yes"}
        } else {intron_gene <- "no"}
        if ("exon" %in% hit_intron_check$class){
          exon <- 1
        }
        if ("start_codon" %in% all_compiled_annotations$class[all_compiled_annotations$id == id]) {
          start_codon <- all_compiled_annotations$start[all_compiled_annotations$id == id & all_compiled_annotations$class == "start_codon"][1]
          stop_codon <- all_compiled_annotations$start[all_compiled_annotations$id == id & all_compiled_annotations$class == "stop_codon"][1]
          # what region of the gene is it in if it is not in the exon? have to account for the instances the gene is on the negative or the positive strand here and if the annotation file even documents exon or intron
          if (start_codon < stop_codon) {
            if (end < start_codon) {upstream_utr <- 1} else if (start > stop_codon) {downstram_utr <- 1} else if ("exon" %in% all_compiled_annotations$class[all_compiled_annotations$id == id]) {intron <- 1}
          } else if (start_codon > stop_codon) {
            if (start > start_codon) {upstream_utr <- 1} else if (end < stop_codon) {downstram_utr <- 1} else if ("exon" %in% all_compiled_annotations$class[all_compiled_annotations$id == id]) {intron <- 1}
          }
        }
      }
    }

    if (i == 1){
      region_metaData <- data.frame("region" = region,
                                    "gene_id" = gene_id,
                                    "intron_gene" = intron_gene,
                                    "gene_bool" = gene_bool,
                                    "Promoter" = Promoter,
                                    "Promoter_id" = Promoter_id,
                                    "exon" = exon,
                                    "intron" = intron,
                                    "CpG_shelf" = CpG_shelf,
                                    "CpG_shore" = CpG_shore,
                                    "CpG_island" = CpG_island
                                    )
      i = 2
    } else {
      region_metaData <- rbind(region_metaData,
                               data.frame("region" = region,
                                          "gene_id" = gene_id,
                                          "intron_gene" = intron_gene,
                                          "gene_bool" = gene_bool,
                                          "Promoter" = Promoter,
                                          "Promoter_id" = Promoter_id,
                                          "exon" = exon,
                                          "intron" = intron,
                                          "CpG_shelf" = CpG_shelf,
                                          "CpG_shore" = CpG_shore,
                                          "CpG_island" = CpG_island
                               ))
    }
  }
  return(region_metaData)
}
