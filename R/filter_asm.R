#'formatting for image
#'
#' @param asm_res_oi the .asm file that you'd like to filter
#' @return returns the snp cpg combinations with allele-specific methylation
#' @export filter_asm

filter_asm <- function(asm_res_oi){
  asm_snpcpg <- unique(asm_res_oi$Chr_SNP_C[asm_res_oi$ASM == TRUE])
  asm_res_oi <- asm_res_oi[asm_res_oi$Chr_SNP_C %in% asm_snpcpg,]
  return(asm_res_oi)
}
