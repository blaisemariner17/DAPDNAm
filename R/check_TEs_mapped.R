#'generate region metadata from regions of interest given provided annotation files
#'
#' @param bam bam filename
#' @param path_to_bams directory path
#' @param DogOverview dog_id metaData
#' @param chrs chromosomes vector
#' @param te_oi_locations locations of interest
#' @return check a bam file for reads overlapping ones of interest
#' @export check_TEs_mapped

check_TEs_mapped <- function(bam, DogOverview=DogOverview, chrs=chrs, te_oi_locations=te_oi_locations, path_to_bams = "/scratch/vsohrab/dap_map_canfam4/nvidia_fq2bam_mapping/bams/"){
  bf <- Rsamtools::BamFile(paste0(path_to_bams, bam))
  dog_id <- gsub(".canfam4.bam","",bam)
  age <- DogOverview$Estimated_Age_Years_at_HLES[DogOverview$dog_id == dog_id]

  fl <- paste0(path_to_bams, bam)

  ## filter to a single file
  param <- Rsamtools::ScanBamParam(
    flag=scanBamFlag(isUnmappedQuery=FALSE),
    what="seq")
  filter <- FilterRules(list(MinWidth = function(x) width(x$seq) > 35))
  dest <- Rsamtools::filterBam(fl, tempfile(), param=param, filter=filter)
  res3 <- Rsamtools::scanBam(dest, param=ScanBamParam(what=c("flag", "qname","rname","pos","qwidth","mapq")))[[1]]

  scan_bam_df <- (data.frame(res3["qname"], res3["flag"], res3["rname"], res3["pos"], res3["qwidth"],res3["mapq"] ))
  #get rid of the PCR or optical duplicates
  scan_bam_df <- scan_bam_df[!scan_bam_df$flag %in% c(1171,1123,1187,1107,1161,1145,1185,1105,1097,1201,1121,1169,1137,1089,1153,1209),]

  scan_bam_df <- scan_bam_df[scan_bam_df$rname %in% chrs,]
  scan_bam_df <- scan_bam_df[! is.na(scan_bam_df$pos),]
  scan_bam_df$start <- scan_bam_df$pos
  scan_bam_df$stop <- scan_bam_df$pos + scan_bam_df$qwidth

  # all_mapped_reads <- nrow(scan_bam_df)
  scan_bam_df<-scan_bam_df[scan_bam_df$mapq >= 1,] # limit the read multimapping a little bit https://sequencing.qcfail.com/articles/mapq-values-are-really-useful-but-their-implementation-is-a-mess/
  mapped_reads <- nrow(scan_bam_df)

  i = 1
  for (chr in chrs){
    # print(chr)
    scan_bam_df_oi <- scan_bam_df[scan_bam_df$rname == paste(chr),]
    te_oi_locations_oi <- te_oi_locations[te_oi_locations$chr == paste(chr),]

    rangesA <- IRanges::IRanges(te_oi_locations_oi$start, te_oi_locations_oi$stop)
    rangesB <- IRanges::IRanges(scan_bam_df_oi$start, scan_bam_df_oi$stop)

    #which regionsB overlap w no regionA regions
    ov <- GenomicRanges::countOverlaps(rangesA, rangesB, type = 'within')

    if (i == 1){
      te_oi_res <- cbind(te_oi_locations_oi, ov)
      i=2
    } else {
      te_oi_res<- rbind(te_oi_res, cbind(te_oi_locations_oi, ov))
    }
  }

  # res_df <- data.frame("dog_id" = dog_id, "Age" = age, "mapped_reads" = mapped_reads, "TE_mapped_reads" = res)

  te_oi_res$dog_id <- dog_id
  te_oi_res$Age <- age
  te_oi_res$mapped_reads <- mapped_reads
  te_oi_res$width <- te_oi_res$stop - te_oi_res$start
  te_oi_res$fpkm <- ((te_oi_res$ov) / (te_oi_res$width * te_oi_res$mapped_reads))*1e10

  return(te_oi_res)
}
