#'generate region, methylation, and coverage data of a bsseq object. can be run in parallel with a list of bsseq obects (e.g. by chromosomes)
#'
#' @param dap_chr_oi the bsseq object to be passed into the function-- can be done in parallel as a list with mclapply
#' @param promoter_annotation promoter annotation of interest. data frame that contains rownames and the columns c("seqnames", "start", "end", "class", "id")
#' @return Function returns methylation and coverage information for the promoter regions passed in
#' @export promoter_coverage_methylation

promoter_coverage_methylation <- function(dap_chr_oi,
                                          promoter_annotation,
                                          at_least_X_coverage_of_a_region = 5,
                                          in_at_least_X_samples = 200){

  chr_oi <- unique(GenomeInfoDb::seqnames(dap_chr_oi))
  regions <- promoter_annotation[promoter_annotation$seqnames == chr_oi,]

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

  return_ <- list(coverage_regions_oi, methylation_regions_oi)
  names(return_)=c("coverage", "methylation")
  return(return_)
}
