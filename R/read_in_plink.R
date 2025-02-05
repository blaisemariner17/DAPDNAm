#'generate plink_in file
#'
#' @param chr chromosome
#' @return Function reads in plink
#' @export read_in_plink_path

read_in_plink_path <- function(chr
) {
  plink_data <- genio::read_plink(paste0("/scratch/bmarine2/dap_rrbs/1_meQTL/meQTL-241101-GEMMA/gemma_out/filtered_genotypes_chr",chr,".bed"))
  plink_in <- plink_data$X
  snp_bims <- data.frame(plink_data$bim)
  snp_bims$chr_snp <- paste0(chr, "_", snp_bims$pos)
  rownames(plink_in) <- snp_bims$chr_snp
  plink_in <- t(plink_in)
  plink_in <- (plink_in[,colSums(plink_in == 2) < .9*nrow(plink_in) & colSums(plink_in == 1) < .9*nrow(plink_in) & colSums(plink_in == 0) < .9*nrow(plink_in)])

  colnames(plink_in) <- gsub(paste0(chr,"_"), "", colnames(plink_in))

  if(chr == "X"){
    colnames(plink_in) <- gsub(paste0("39_"), "", colnames(plink_in))
  }

  rownames(plink_in) <- gsub(".canfam4", "", rownames(plink_in))
  rm(plink_data)
  # Combine the results into a dataframe
  plink_in<- data.frame(t(unique(t(plink_in))))
  colnames(plink_in) <- gsub("X", "", colnames(plink_in))

  if(chr=="X"){
    plink_in<- plink_in[rownames(plink_in) %in% metaData$dog_id,]

    plink_data <- genio::read_plink(paste0("/scratch/bmarine2/dap_rrbs/1_meQTL/meQTL-241101-GEMMA/gemma_out/filtered_genotypes_notPAR_chrX.bed"))
    plink_in2 <- plink_data$X
    snp_bims <- data.frame(plink_data$bim)
    snp_bims$chr_snp <- paste0(chr, "_", snp_bims$pos)
    rownames(plink_in2) <- snp_bims$chr_snp
    plink_in2 <- t(plink_in2)
    plink_in2 <- (plink_in2[,colSums(plink_in2 == 2) < .9*nrow(plink_in2) & colSums(plink_in2 == 1) < .9*nrow(plink_in2) & colSums(plink_in2 == 0) < .9*nrow(plink_in2)])

    colnames(plink_in2) <- gsub(paste0(chr,"_"), "", colnames(plink_in2))

    if(chr == "X"){
      colnames(plink_in2) <- gsub(paste0("39_"), "", colnames(plink_in2))
    }

    rownames(plink_in2) <- gsub(".canfam4", "", rownames(plink_in2))
    rm(plink_data)
    # Combine the results into a dataframe
    plink_in2<- data.frame(t(unique(t(plink_in2))))
    colnames(plink_in2) <- gsub("X_", "", colnames(plink_in2))
    colnames(plink_in2) <- gsub("X", "", colnames(plink_in2))

    plink_in2<- plink_in2[rownames(plink_in2) %in% metaData$dog_id,]

    plink_in <- plink_in[rownames(plink_in) %in% rownames(plink_in2),]
    plink_in2 <- plink_in2[rownames(plink_in2) %in% rownames(plink_in),]

    plink_in <- plink_in[order(rownames(plink_in)),]
    plink_in2 <- plink_in2[order(rownames(plink_in2)),]

    test <- cbind(plink_in, plink_in2)
    plink_in <- test
    rm(test)

  }
  return(plink_in)
}
