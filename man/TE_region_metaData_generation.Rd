% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TE_region_metaData_generation.R
\name{TE_region_metaData_generation}
\alias{TE_region_metaData_generation}
\title{generate transposable element region metadata from regions of interest given provided annotation files}
\usage{
TE_region_metaData_generation(
  regions,
  transposons_annotation,
  genome_gene_annotation,
  sw_score_cutoff = 225,
  expanded_search_if_nec = 50,
  return_classes = FALSE
)
}
\arguments{
\item{regions}{regions of interest}

\item{transposons_annotation}{TE annotation file}

\item{genome_gene_annotation}{gene annotation to find the nearest gene}

\item{sw_score_cutoff}{the minimum sw score that will be queried here}

\item{expanded_search_if_nec}{expanding the search if necessary for region overlaps}

\item{return_classes}{if TRUE, returns if the region overlaps with transposon classes (LINE1, LINE2, SINE_tRNA, etc.)}
}
\value{
Function returns region metadata of interest
}
\description{
generate transposable element region metadata from regions of interest given provided annotation files
}
