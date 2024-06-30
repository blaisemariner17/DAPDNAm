#'generate region metadata from regions of interest given provided annotation files
#'
#' @param bam bam filename
#' @param path_to_bams directory path
#' @param DogOverview dog_id metaData
#' @param chrs chromosomes vector
#' @param locations_oi locations of interest
#' @return check a bam file for reads overlapping ones of interest
#' @export check_TEs_mapped

check_TEs_mapped <- function(bam, DogOverview=DogOverview, chrs=chrs, te_oi_locations=te_oi_locations, path_to_bams = "/scratch/vsohrab/dap_map_canfam4/nvidia_fq2bam_mapping/bams/"){
  bf <- Rsamtools::BamFile(paste0(path_to_bams, bam))
  dog_id <- gsub(".canfam4.bam","",bam)
  age <- DogOverview$Estimated_Age_Years_at_HLES[DogOverview$dog_id == dog_id]

  fl <- paste0(path_to_bams, bam)

  ## filter to a single file
  param <- ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
  filter <- FilterRules(list(MinWidth = function(x) width(x$seq) > 35))
  dest <- filterBam(fl, tempfile(), param=param, filter=filter)
  res3 <- scanBam(dest, param=ScanBamParam(what=c("qname","rname","pos","qwidth","mapq")))[[1]]

  scan_bam_df <- (data.frame(res3["qname"], res3["rname"], res3["pos"], res3["qwidth"],res3["mapq"] ))

  scan_bam_df <- scan_bam_df[scan_bam_df$rname %in% chrs,]
  scan_bam_df <- scan_bam_df[! is.na(scan_bam_df$pos),]
  scan_bam_df$start <- scan_bam_df$pos
  scan_bam_df$stop <- scan_bam_df$pos + scan_bam_df$qwidth

  all_mapped_reads <- nrow(scan_bam_df)
  scan_bam_df<-scan_bam_df[scan_bam_df$mapq >=20,] #99% chance of unique
  uniquely_mapped_reads <- nrow(scan_bam_df)

  res <- 0
  for (chr in chrs){
    scan_bam_df_oi <- scan_bam_df[scan_bam_df$rname == paste(chr),]
    te_oi_locations_oi <- te_oi_locations[te_oi_locations$chr == paste(chr),]

    rangesA <- IRanges::IRanges(te_oi_locations_oi$start, te_oi_locations_oi$stop)
    rangesB <- IRanges::IRanges(scan_bam_df_oi$start, scan_bam_df_oi$stop)

    #which regionsB overlap w no regionA regions
    ov <- GenomicRanges::countOverlaps(rangesB, rangesA, type = 'any')>=1
    res <- res + sum(ov)
  }


  # print(paste("Age:", age, ". TE overlaps:", res))

  res_df <- data.frame("dog_id" = dog_id, "Age" = age, "all_mapped_reads" = all_mapped_reads, "uniquely_mapped_reads" = uniquely_mapped_reads, "TE_mapped_reads" = res)

  return(res_df)
}
