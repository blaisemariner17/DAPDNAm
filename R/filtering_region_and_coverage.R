#'generate region, methylation, and coverage data of a bsseq object. can be run in parallel with a list of bsseq obects (e.g. by chromosomes)
#'
#' @param dap_chr_oi the bsseq object to be passed into the function-- can be done in parallel as a list with mclapply
#' @param metaData sample metaData for colnames of the dap_list
#' @param maxgap keyword arg for region definition for regionFinder3
#' @param at_least_n_sites_per_region lower bound (inclusive) of CpG sites for a region to be considered in this analysis (default = 5)
#' @param at_least_X_coverage_of_a_region used in concert with in_at_least_X_samples (default 5x coverage in 200 samples)
#' @param in_at_least_X_samples used in concert with at_least_X_coverage_of_a_region (default 5x coverage in 200 samples)
#' @param upper_bound omit regions with greater than X mean percent methylation (default = 0.9, inclusive)
#' @param lower_bound omit regions with less than X mean percent methylation (default = 0.1, inclusive)
#' @param regions_manual_input provide own regions if you want, optional
#' @return Function returns region, methylation, and coverage information for the regions sought by regionFinder3
#' @export filtering_region_and_coverage

filtering_region_and_coverage <- function(dap_chr_oi,
                                          metaData = NULL,
                                          maxgap = 50,
                                          at_least_n_sites_per_region = 5,
                                          at_least_X_coverage_of_a_region = 10,
                                          in_at_least_X_samples = 200,
                                          lower_bound = 0.1,
                                          upper_bound = 0.9,
                                          regions_manual_input = NULL){

  if (is.null(regions_manual_input)){
    p1_metaData <- metaData[metaData$first_rrbs == "yes" & grepl("precision", metaData$Cohort),]
    #get the regions
    dap_chr_oi_p1 <- dap_chr_oi[,sampleNames(dap_chr_oi) %in% p1_metaData$lid_pid]
    regions <- bsseq:::regionFinder3(x = as.integer(rep(1,
                                                        length(dap_chr_oi_p1))),
                                     chr = as.character(GenomeInfoDb::seqnames(dap_chr_oi_p1)),
                                     positions = BiocGenerics::start(dap_chr_oi_p1), maxGap = maxgap,
                                     verbose = FALSE)[["up"]]
    rownames(regions) <- paste(regions$chr, regions$start, regions$end, sep = "_")
    regions <- regions[regions$n >= at_least_n_sites_per_region,]
    regions$length <- regions$end - regions$start + 1

  } else {
    regions <-  regions_manual_input[startsWith(rownames(regions_manual_input), prefix = paste0(unique(as.character(GenomeInfoDb::seqnames(dap_chr_oi))),"_")),]
  }

  coverage_regions_oi <- bsseq::getCoverage(dap_chr_oi, type = "Cov", regions = regions,
                                            what="perRegionTotal", withDimnames = TRUE)

  rownames(coverage_regions_oi) <- rownames(regions)

  if (is.null(regions_manual_input)){
    coverage_regions_oi <- coverage_regions_oi[rowSums(coverage_regions_oi >= at_least_X_coverage_of_a_region) >= in_at_least_X_samples, ]
  }

  methylation_regions_oi <- bsseq::getCoverage(dap_chr_oi, type = "M", regions = regions,
                                               what="perRegionTotal", withDimnames = TRUE)

  rownames(methylation_regions_oi) <- rownames(regions)

  methylation_regions_oi <- methylation_regions_oi[rownames(methylation_regions_oi) %in% rownames(coverage_regions_oi),]

  if (is.null(regions_manual_input)){
    regions_oi <- regions[rownames(regions) %in% rownames(methylation_regions_oi),]

    perc_meth <- methylation_regions_oi / coverage_regions_oi
    perc_meth <- perc_meth[matrixStats::rowMedians(perc_meth, na.rm = TRUE) >= lower_bound & matrixStats::rowMedians(perc_meth, na.rm = TRUE) <= upper_bound,]
    regions_oi <- regions[rownames(regions) %in% rownames(perc_meth),]
    coverage_regions_oi <- coverage_regions_oi[rownames(coverage_regions_oi) %in% rownames(perc_meth),]
    methylation_regions_oi <- methylation_regions_oi[rownames(methylation_regions_oi) %in% rownames(perc_meth),]
  } else {
    regions_oi <- regions
  }

  return_ <- list(regions_oi, coverage_regions_oi, methylation_regions_oi)
  names(return_)=c("regions", "coverage", "methylation")
  return(return_)
}
