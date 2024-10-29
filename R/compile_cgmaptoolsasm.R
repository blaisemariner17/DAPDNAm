#'compile the cgmaptools asm files
#'
#' @param chr chromosome oi (e.g., 'chr1')
#' @param path path to asm files
#' @return returns the compiled asm files for the chromosome of interest
#' @export compile_cgmaptoolsasm

compile_cgmaptoolsasm <- function(chr, path = "/scratch/bmarine2/meQTL-241004/traditional_asm_processing") {
  print(chr)
  NUMIDS <- length(list.files(paste(path), pattern = paste0(chr, ".asm")))

  i=1
  for(file_oi in list.files(paste(path), pattern = paste0(chr, ".asm"))){

    if(file.size(paste0(path, file_oi)) == 0L){next}

    lid_pid<-gsub(paste0(".", chr, ".asm"), "" , file_oi)
    # print(lid_pid)

    asm_oi <- read.table(paste0(path, file_oi), header = T)
    asm_oi$lid_pid <- lid_pid

    if(i==1){asm_chr<-asm_oi;i=2}else{asm_chr<-rbind(asm_chr,asm_oi)}
  }
  if(!dir.exists("asms_compiled")){dir.create("asms_compiled")}
  write.table(asm_chr, paste0("asms_compiled/", chr,"_ASM.ASM"))

  #remove the cpgs that are in less than 10% of individuals (Fan et al., 2019)
  asm_chr <- asm_chr[asm_chr$C_Pos %in% names(table(asm_chr$C_Pos)[table(asm_chr$C_Pos) > .1*NUMIDS]),]

  #remove the snp cpg instances found with >5 total reads (Fan et al., 2019)
  asm_chr$Allele1_reads <- as.numeric(gsub("-.*","", asm_chr$Allele1_linked_C)) + as.numeric(gsub(".*-","", asm_chr$Allele1_linked_C))
  asm_chr$Allele2_reads <- as.numeric(gsub("-.*","", asm_chr$Allele2_linked_C)) + as.numeric(gsub(".*-","", asm_chr$Allele2_linked_C))
  asm_chr$total_reads <- asm_chr$Allele1_reads + asm_chr$Allele2_reads
  asm_chr <- asm_chr[asm_chr$total_reads >=5,]

  #remove the perpetually methylated snp cpgs >90% or <10% (Fan et al., 2019)
  asm_chr$Allele1_meth <- as.numeric(gsub("-.*","", asm_chr$Allele1_linked_C))
  asm_chr$Allele2_meth <- as.numeric(gsub("-.*","", asm_chr$Allele2_linked_C))

  asm_chr$Allele1_pmeth <- asm_chr$Allele1_meth / asm_chr$Allele1_reads
  asm_chr$Allele2_pmeth <- asm_chr$Allele2_meth / asm_chr$Allele2_reads

  asm_chr <- asm_chr[asm_chr$Allele1_pmeth < 0.9 | asm_chr$Allele1_pmeth > 0.1 |
                       asm_chr$Allele2_pmeth < 0.9 | asm_chr$Allele2 > 0.1,]

  # require at least 5% mean allele frequency up to here (Fan et al., 2019)
  asm_chr <- asm_chr[asm_chr$SNP_Pos %in% names(table(asm_chr$SNP_Pos)[table(asm_chr$SNP_Pos) > .05*NUMIDS]),]

  #label the geno of the snp
  # and their "geno" variable is also infuriating but 0 means the reference allele matches both allele 1 and allele 2, 2 means neither match, and 1 means else....
  asm_chr$geno <- 1
  asm_chr$geno[asm_chr$Ref == asm_chr$Allele1 & asm_chr$Allele1 == asm_chr$Allele2] <- 0
  asm_chr$geno[asm_chr$Ref != asm_chr$Allele1 & asm_chr$Allele1 != asm_chr$Allele2] <- 2

  asm_chr$chr <- chr
  asm_chr$Chr_SNP_C <- paste(asm_chr$chr, asm_chr$SNP_Pos, asm_chr$C_Pos, sep = "_")

  #total meth reads
  asm_chr$total_meth <- asm_chr$Allele1_meth + asm_chr$Allele1_meth

  return(asm_chr)
}
