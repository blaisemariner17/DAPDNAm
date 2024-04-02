#' clock simulations
#'
#' @param alphs_oi alphas to use for glmnet
#' @param metaData sample metaData
#' @param region_metaData region metaData
#' @param perc_meth imputed percent methylation matrix
#' @param region_metaData_te TE region metaData
#' @param perc_meth_te imputed TE percent methylation matrix
#' @param region_metaData_prmtr promoter region metaData
#' @param perc_meth_prmtr imputed promoter percent methylation matrix
#' @param region_metaData_cpgisl promoter region metaData
#' @param perc_meth_cpgisl imputed promoter percent methylation matrix
#' @return Function returns a list of the simulated clock results
#' @export clock_simulate

clock_simulate <- function(alph,
                           metaData,
                           region_metaData,
                           perc_meth,
                           region_metaData_te,
                           perc_meth_te,
                           region_metaData_prmtr,
                           perc_meth_prmtr,
                           region_metaData_cpgisl,
                           perc_meth_cpgisl
) {
  # print(i)
  sampling <- sample(rownames(perc_meth), as.integer(nrow(perc_meth_te)))
  perc_meth_rand <- perc_meth[rownames(perc_meth) %in% sampling,]
  region_metaData_rand <- region_metaData[region_metaData$region %in% sampling,]

  sampling <- sample(unique(metaData$dog_id), 100)

  list_res <- DAPDNAm::build_clock(alph,
                                   metaData,
                                   region_metaData,
                                   perc_meth,
                                   test_samples = colnames(perc_meth) %in% metaData$lid_pid[metaData$dog_id %in% sampling]
  )

  list_res_TE <- DAPDNAm::build_clock(alph,
                                      metaData,
                                      region_metaData_te,
                                      perc_meth_te,
                                      test_samples = colnames(perc_meth_te) %in% metaData$lid_pid[metaData$dog_id %in% sampling]
  )

  list_res_PR <- DAPDNAm::build_clock(alph,
                                      metaData,
                                      region_metaData_prmtr,
                                      perc_meth_prmtr,
                                      test_samples = colnames(perc_meth_prmtr) %in% metaData$lid_pid[metaData$dog_id %in% sampling]
  )

  list_res_CpGIs <- DAPDNAm::build_clock(alph,
                                         metaData,
                                         region_metaData_cpgisl,
                                         perc_meth_cpgisl,
                                         test_samples = colnames(perc_meth_cpgisl) %in% metaData$lid_pid[metaData$dog_id %in% sampling]
  )

  list_res_rand <- DAPDNAm::build_clock(alph,
                                        metaData,
                                        region_metaData_rand,
                                        perc_meth_rand,
                                        test_samples = colnames(perc_meth_rand) %in% metaData$lid_pid[metaData$dog_id %in% sampling]
  )

  list_res[['metaData']]$regions <- "All"
  list_res_TE[['metaData']]$regions <- "TE"
  list_res_PR[['metaData']]$regions <- "Promoter"
  list_res_CpGIs[['metaData']]$regions <- "CpG_island"
  list_res_rand[['metaData']]$regions <- "Random"

  list_res[['region_metaData']]$regions <- "All"
  list_res_TE[['region_metaData']]$regions <- "TE"
  list_res_PR[['region_metaData']]$regions <- "Promoter"
  list_res_CpGIs[['region_metaData']]$regions <- "CpG_island"
  list_res_rand[['region_metaData']]$regions <- "Random"

  list_res[['clock_results']]$regions <- "All"
  list_res_TE[['clock_results']]$regions <- "TE"
  list_res_PR[['clock_results']]$regions <- "Promoter"
  list_res_CpGIs[['clock_results']]$regions <- "CpG_island"
  list_res_rand[['clock_results']]$regions <- "Random"

  all_res_list <- list()
  # all_res_list[["metaData"]] <- rbind(list_res[['metaData']], list_res_TE[['metaData']],list_res_PR[['metaData']], list_res_CpGIs[['metaData']],list_res_rand[['metaData']])
  # all_res_list[["region_metaData"]] <- rbind(list_res[['region_metaData']], list_res_TE[['region_metaData']],list_res_PR[['region_metaData']], list_res_CpGIs[['region_metaData']],list_res_rand[['region_metaData']])
  all_res_list[["clock_results"]] <- rbind(list_res[['clock_results']], list_res_TE[['clock_results']],list_res_PR[['clock_results']], list_res_CpGIs[['clock_results']],list_res_rand[['clock_results']])

  return(all_res_list)
}
