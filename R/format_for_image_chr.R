#'formatting for image by chr
#'
#' @param asm_chr the .asm file that you'd like to filter
#' @return returns the image formated requirements
#' @export format_for_image_chr

VectorIntersect <- function(v,z) {
  unlist(lapply(unique(v[v%in%z]), function(x) rep(x,min(sum(v==x),sum(z==x)))))
}
is.contained <- function(v,z) {length(VectorIntersect(v,z))==length(v)}

format_for_image_chr <- function(asm_chr){

  chr <- unique(asm_chr$Chr)

  asm_chr <- asm_chr[order(asm_chr$lid_pid),]
  SNP_CPG <- unique(gsub("_L.*", "", asm_chr$Chr_SNP_C))

  asm_chr$lid_pid <- gsub(paste0("_", chr), "", asm_chr$lid_pid)
  lid_pids <- unique(asm_chr$lid_pid)

  for (lid_pid in lid_pids){
    asm_chr_lidpid <- asm_chr[asm_chr$lid_pid == lid_pid,]

    #add all the snps not in the lidpid
    if (!is.contained(SNP_CPG,asm_chr_lidpid$Chr_SNP_C)){
      not_in <- SNP_CPG[!SNP_CPG %in% asm_chr_lidpid$Chr_SNP_C]

      remove_chr <- gsub(paste0(chr, "_"), "", not_in)
      SNP_Pos <- gsub("_.*","",remove_chr)
      C_Pos <- gsub(".*_","",remove_chr)

      rbinding <- data.frame("Chr_SNP_C" = not_in,
                             "SNP_Pos" = SNP_Pos,
                             "C_Pos" = C_Pos,
                             "lid_pid" = rep(lid_pid, length(not_in))
      )

      for (col in colnames(asm_chr_lidpid)[!colnames(asm_chr_lidpid) %in% colnames(rbinding)]){
        rbinding[,c(paste(col))] <- 0
      }
      asm_chr_lidpid <- rbind(asm_chr_lidpid,rbinding)
    }
    asm_chr_lidpid <- asm_chr_lidpid[order(asm_chr_lidpid$Chr_SNP_C),]

    if(lid_pid == lid_pids[1]){
      genotype <- asm_chr_lidpid[,c("Chr_SNP_C", "geno")]
      r <- asm_chr_lidpid[,c("Chr_SNP_C", "total_reads")]
      r1 <- asm_chr_lidpid[,c("Chr_SNP_C", "Allele1_reads")]
      r2 <- asm_chr_lidpid[,c("Chr_SNP_C", "Allele2_reads")]
      y <- asm_chr_lidpid[,c("Chr_SNP_C", "total_meth")]
      y1 <- asm_chr_lidpid[,c("Chr_SNP_C", "Allele1_meth")]
      y2 <- asm_chr_lidpid[,c("Chr_SNP_C", "Allele2_meth")]

      colnames(genotype) <- c("Chr_SNP_C", paste(lid_pid))
      colnames(r) <- c("Chr_SNP_C", paste(lid_pid))
      colnames(r1) <- c("Chr_SNP_C", paste(lid_pid))
      colnames(r2) <- c("Chr_SNP_C", paste(lid_pid))
      colnames(y) <- c("Chr_SNP_C", paste(lid_pid))
      colnames(y1) <- c("Chr_SNP_C", paste(lid_pid))
      colnames(r2) <- c("Chr_SNP_C", paste(lid_pid))
      colnames(y2) <- c("Chr_SNP_C", paste(lid_pid))
    } else {
      genotype[,c(paste(lid_pid))] <- asm_chr_lidpid$geno
      r[,c(paste(lid_pid))] <- asm_chr_lidpid$total_reads
      r1[,c(paste(lid_pid))] <- asm_chr_lidpid$Allele1_reads
      r2[,c(paste(lid_pid))] <- asm_chr_lidpid$Allele2_reads
      y[,c(paste(lid_pid))] <- asm_chr_lidpid$total_meth
      y1[,c(paste(lid_pid))] <- asm_chr_lidpid$Allele1_meth
      y2[,c(paste(lid_pid))] <- asm_chr_lidpid$Allele2_meth

      if(lid_pid == lid_pids[2]){
        rownames(genotype) <- genotype$Chr_SNP_C
        rownames(r)<- r$Chr_SNP_C
        rownames(r1)<- r1$Chr_SNP_C
        rownames(r2)<- r2$Chr_SNP_C
        rownames(y)<- y$Chr_SNP_C
        rownames(y1)<- y1$Chr_SNP_C
        rownames(y2)<- y2$Chr_SNP_C

        genotype <- genotype[,colnames(genotype) != "Chr_SNP_C"]
        r <- r[,colnames(r) != "Chr_SNP_C"]
        y <- y[,colnames(y) != "Chr_SNP_C"]
        r1 <- r1[,colnames(r1) != "Chr_SNP_C"]
        r2 <- r2[,colnames(r2) != "Chr_SNP_C"]
        y1 <- y1[,colnames(y1) != "Chr_SNP_C"]
        y2 <- y2[,colnames(y2) != "Chr_SNP_C"]
      }

    }
  }

  format_for_image_chr <- list()
  format_for_image_chr[['geno']] <- (genotype)
  format_for_image_chr[['r']] <- t(r)
  format_for_image_chr[['y']] <- t(y)
  format_for_image_chr[['r1']] <- t(r1)
  format_for_image_chr[['r2']] <- t(r2)
  format_for_image_chr[['y1']] <- t(y1)
  format_for_image_chr[['y2']] <- t(y2)
  format_for_image_chr[['lid_pids']] <- lid_pids

  return(format_for_image_chr)
}
