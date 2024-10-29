#'
#' @param format_for_image_chr the image formatted data
#' @return returns image results
#' @export IMAGE_fnct

IMAGE_fnct <- function(format_for_image_chr){
  genotype <- format_for_image_chr[["geno"]]
  lid_pids <- format_for_image_chr[['lid_pids']]

  Kin <- readRDS("../BASELINE-PRECISION-240820/metaData_samples/240821-baseline_GRM.rds")
  Kin <- Kin[rownames(Kin) %in% lid_pids, colnames(Kin) %in% lid_pids]

  geno<-list()
  geno[[1]]<-matrix(0,ncol=ncol(genotype),nrow=nrow(genotype))
  geno[[2]]<-matrix(0,ncol=ncol(genotype),nrow=nrow(genotype))
  for(i in 1:ncol(genotype)){
    for(j in 1:nrow(genotype))
    {
      if(is.na(genotype[j,i]))
      {
        geno[[2]][j,i]=0/0
        geno[[1]][j,i]=0/0
      }else if(genotype[j,i]==1){
        geno[[2]][j,i]=1
      }else if(genotype[j,i]==2){
        geno[[2]][j,i]=1
        geno[[1]][j,i]=1
      }
    }
  }
  names(geno) <- c('hap1', 'hap2')

  data <- list()
  data[[1]]<-format_for_image_chr[['r']]
  data[[2]]<-format_for_image_chr[['y']]
  data[[3]]<-format_for_image_chr[['r1']]
  data[[4]]<-format_for_image_chr[['r2']]
  data[[5]]<-format_for_image_chr[['y1']]
  data[[6]]<-format_for_image_chr[['y2']]
  names(data) <- c('r', 'y', 'r1', 'r2', 'y1', 'y2')

  res <- IMAGE::image(geno,
                      data,
                      Kin,
                      # Covariates = cov_573_int,
                      numCore=1
  )

  tested <- res$Loc
  res$Chr_SNP_C <- colnames(data[['y2']])[tested]

  return(res)
}
