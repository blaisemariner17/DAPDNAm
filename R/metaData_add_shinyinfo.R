#'build sample metadata
#'
#' @param metaData alpha to use for glmnet
#' @param shiny_metaData alpha to use for glmnet
#' @param cols_oi what cols to cbind to the metaData
#' @return Function returns the metadata with columns from shiny
#' @export shiny_meta_add

shiny_meta_add <- function(metaData, shiny_metaData, cols_oi = c("rdid", "prep_date")) {
  metaData_oi <- metaData[metaData$sid %in% shiny_metadata$dap_sample_id,]
  shiny_oi <- shiny_metaData[shiny_metaData$dap_sample_id %in% metaData_oi$sid,]
  metaData_oi <- metaData_oi[order(metaData_oi$sid),]
  shiny_oi <- shiny_oi[order(shiny_oi$dap_sample_id),]
  if (all(shiny_oi$dap_sample_id == metaData_oi$sid)) {
    metaData_oi <- cbind(metaData_oi, shiny_oi[,cols_oi])
    return(metaData_oi)
  } else {warning("check inputs")}
}
