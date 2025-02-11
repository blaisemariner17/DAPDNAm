% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/promoter_coverage_methylation.R
\name{promoter_coverage_methylation}
\alias{promoter_coverage_methylation}
\title{generate region, methylation, and coverage data of a bsseq object. can be run in parallel with a list of bsseq obects (e.g. by chromosomes)}
\usage{
promoter_coverage_methylation(
  dap_chr_oi,
  promoter_annotation,
  at_least_X_coverage_of_a_region = 5,
  in_at_least_X_samples = 200
)
}
\arguments{
\item{dap_chr_oi}{the bsseq object to be passed into the function-- can be done in parallel as a list with mclapply}

\item{promoter_annotation}{promoter annotation of interest. data frame that contains rownames and the columns c("seqnames", "start", "end", "class", "id")}
}
\value{
Function returns methylation and coverage information for the promoter regions passed in
}
\description{
generate region, methylation, and coverage data of a bsseq object. can be run in parallel with a list of bsseq obects (e.g. by chromosomes)
}
