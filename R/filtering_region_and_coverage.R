#'generate region, methylation, and coverage data of a bsseq object. can be run in parallel with a list of bsseq obects (e.g. by chromosomes)
#'
#' @param dap_chr_oi the bsseq object to be passed into the function-- can be done in parallel as a list with mclapply
#' @param upper_bound omit regions with greater than X mean percent methylation
#' @param lower_bound omit regions with less than X mean percent methylation
#' @return Function returns region, methylation, and coverage information for the regions sought by regionFinder3
#' @export filtering_region_and_coverage

filtering_region_and_coverage <- function(dap_chr_oi,
                                          lower_bound = 0.1,
                                          upper_bound = 0.9){
  regions <- bsseq:::regionFinder3(x = as.integer(rep(1,
                                                      length(dap_chr_oi))),
                                   chr = as.character(GenomeInfoDb::seqnames(dap_chr_oi)),
                                   positions = BiocGenerics::start(dap_chr_oi), maxGap = maxgap,
                                   verbose = FALSE)[["up"]]

  rownames(regions) <- paste(regions$chr, regions$start, regions$end, sep = "_")

  regions <- regions[regions$n >= 5,]
  regions$length <- regions$end - regions$start

  coverage_regions_oi <- bsseq::getCoverage(dap_chr_oi, type = "Cov", regions = regions,
                                            what="perRegionTotal", withDimnames = TRUE)

  rownames(coverage_regions_oi) <- rownames(regions)

  coverage_regions_oi <- coverage_regions_oi[rowSums(coverage_regions_oi >= 5) > 200, ]

  methylation_regions_oi <- bsseq::getCoverage(dap_chr_oi, type = "M", regions = regions,
                                               what="perRegionTotal", withDimnames = TRUE)

  rownames(methylation_regions_oi) <- rownames(regions)

  methylation_regions_oi <- methylation_regions_oi[rownames(methylation_regions_oi) %in% rownames(coverage_regions_oi),]

  regions_oi <- regions[rownames(regions) %in% rownames(methylation_regions_oi),]

  perc_meth <- methylation_regions_oi / coverage_regions_oi
  perc_meth <- perc_meth[rowMeans(perc_meth) >= lower_bound & rowMeans(perc_meth) <= upper_bound,]
  regions_oi <- regions_all_chr[rownames(regions_all_chr) %in% rownames(perc_meth),]
  coverage_regions_oi <- coverage_regions_oi[rownames(coverage_regions_oi) %in% rownames(perc_meth),]
  methylation_regions_oi <- methylation_regions_oi[rownames(methylation_regions_oi) %in% rownames(perc_meth),]

  return_ <- list(regions_oi, coverage_regions_oi, methylation_regions_oi)
  names(return_)=c("regions", "coverage", "methylation")
  return(return_)
}
