#### Blaise Mariner
### this script generates the region metaData such that there is a way to better understand the reigons that are differentially methylated with different demographics
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-02-13
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project

#### GLOBALS ####
# this scrip requires at least 12 cores
cores_ <- 20

#### END GLOBALS ####

library_list <- c(
  "ggplot2",
  "svglite",
  "ggExtra",
  "ggtext",
  "GenomicRanges",
  "GenomicFeatures",
  "tidyverse",
  "bsseq",
  "comethyl",
  "PQLseq",
  "foreach",
  "parallel",
  "Biostrings",
  "patchwork"
)

lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 24, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      axis.text = element_text(family = "sans", size = 18, color="black"),       axis.title = element_text(family = "sans", size = 20, color="black"),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),      legend.text.align =	0,      panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

gtf <- rtracklayer::import('/home/bmarine2/GENOME-ANNOTATION-FILE/UU_Cfam_GSD_1.0_ROSY.refSeq.ensformat.gtf')
gtf_cpgisl <- read.csv("/home/bmarine2/GENOME-ANNOTATION-FILE/cpgIslandExt.txt", header = F, sep = "\t")
bed_file <- "/home/bmarine2/GENOME-ANNOTATION-FILE/SP_13_dense.bed"
chromatin_states_annotation = rtracklayer::import(bed_file, format = "bed")
chromatin_states_annotation$name <- as.character(chromatin_states_annotation$name)
chromatin_states_annotation$name[chromatin_states_annotation$name == '1' | chromatin_states_annotation$name == '2'|
                                   chromatin_states_annotation$name == '3' | chromatin_states_annotation$name == '4'] <- "ChrSt_promoter"
chromatin_states_annotation$name[chromatin_states_annotation$name == '5' | chromatin_states_annotation$name == '6'|
                                   chromatin_states_annotation$name == '7'] <- "ChrSt_enhancer"
chromatin_states_annotation$name[chromatin_states_annotation$name == '8' | chromatin_states_annotation$name == "9"] <- "ChrSt_polycomb"
chromatin_states_annotation$name[chromatin_states_annotation$name == '10' | chromatin_states_annotation$name == "11"|
                                   chromatin_states_annotation$name == "12"] <- "ChrSt_heterochromatin"
chromatin_states_annotation$name[chromatin_states_annotation$name == "13"] <- "ChrSt_quies"

gtf_trans <- read.csv("/home/bmarine2/GENOME-ANNOTATION-FILE/rmsk.txt", header = F, sep = "\t")
colnames(gtf_trans) <- c('bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns'	,
                         'seqnames', 'start', 'end', 'genoLeft', 'strand'	,
                         'repName', 'class', 'repFamily', 'repStart', 'repEnd',	'repLeft',	'id')
gtf_trans <- gtf_trans[gtf_trans$seqnames %in% paste0("chr", c(1:38, "X")),]

# gtf_trans$length<-gtf_trans$end-gtf_trans$start
# for_gg <- gtf_trans[gtf_trans$class == "LINE",]
# for_gg <- for_gg[grepl("L1M", for_gg$repName),]
# ggplot(for_gg,aes(x=repName, y = length))+
#   geom_boxplot()+theme_blaise+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


transposon_annotation <- gtf_trans[gtf_trans$class %in% c("LINE", "SINE", "LTR", "DNA",
                                                          "Simple_repeat","Low_complexity", "rRNA", "tRNA",
                                                          "Satellite","Unknown"
),]

dim(transposon_annotation)

rm(gtf_trans)

# .rs.restartR()
# remove.packages("DAPDNAm")
if (! "DAPDNAm" %in% installed.packages()){devtools::install_github("blaisemariner17/DAPDNAm", force = TRUE)}
library("DAPDNAm")
#get the maxgaps I used from the previous loop
maxgaps_split <- str_split(list.files(path = "getcoverage_data_and_plots/", pattern = "_all_coverage_regions_oi.rds"), fixed("_"))
maxgaps <- c()
for (i in 1:length(maxgaps_split)) {
  maxgap <- maxgaps_split[[i]][2]
  maxgaps <- append(maxgaps, maxgap)
}
maxgaps = "250"
if(dir.exists("metadata_regions")==F){dir.create("metadata_regions")}

for (maxgap in maxgaps){
  print(maxgap)
  regions <- read_rds(paste0("getcoverage_data_and_plots/maxGap_",maxgap,"_all_regions_oi.rds"))

  chrs=paste0("chr",c(1:38,"X"))
  regions_list <- parallel::mclapply(chrs,
                                     function(x, regions){return(regions[regions$chr == x,])},
                                     regions = regions,
                                     mc.cores =18)
  names(regions_list) <- chrs
  region_metaData_TE_list <- parallel::mclapply(regions_list,
                                                DAPDNAm::TE_region_metaData_generation,
                                                transposons_annotation = transposon_annotation,
                                                sw_score_cutoff = 225,
                                                mc.cores = 20)

  for (i in 1:length(region_metaData_TE_list)) {
    if (i == 1){region_metaData_TE <- region_metaData_TE_list[[i]]} else {region_metaData_TE <- rbind(region_metaData_TE, region_metaData_TE_list[[i]])}
  }

  region_metaData_list <- parallel::mclapply(regions_list,
                                             DAPDNAm::region_metaData_generation,
                                             genome_gene_annotation = gtf,
                                             cpgisland_annotation = gtf_cpgisl,
                                             mc.cores = 12)

  for (i in 1:length(region_metaData_list)) {
    if (i == 1){region_metaData <- region_metaData_list[[i]]} else {region_metaData <- rbind(region_metaData, region_metaData_list[[i]])}
  }

  region_metaData_ChromStates_list <- parallel::mclapply(regions_list,
                                                         DAPDNAm::chromatin_state_region_metaData_generation,
                                                         chromatin_states_annotation = chromatin_states_annotation,
                                                         mc.cores = 12)

  for (i in 1:length(region_metaData_ChromStates_list)) {
    if (i == 1){region_metaData_ChromStates <- region_metaData_ChromStates_list[[i]]} else {region_metaData_ChromStates <- rbind(region_metaData_ChromStates, region_metaData_ChromStates_list[[i]])}
  }

  if (all(rownames(region_metaData_TE) == rownames(region_metaData))){region_metaData <- cbind(region_metaData, region_metaData_TE)} else {print("not all the transposable element rownames are the same!!!")}
  if (all(rownames(region_metaData_ChromStates) == rownames(region_metaData))){region_metaData <- cbind(region_metaData, region_metaData_ChromStates)} else {print("not all the chromatin states rownames are the same!!!")}
  # if (all(rownames(human_conserved_region_metaData) == rownames(region_metaData))){region_metaData <- cbind(region_metaData, human_conserved_region_metaData)} else {print("not all the human conserved regions rownames are the same!!!")}

  write_rds(x = region_metaData, file = paste0("metadata_regions/metaData_regions_maxgap_",maxgap,".rds"))

  coverage_all_chr <- read_rds(paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_coverage_regions_oi.rds"))
  methylation_all_chr <- read_rds(paste0("getcoverage_data_and_plots/maxGap_", maxgap, "_all_methylation_regions_oi.rds"))

  plots <- DAPDNAm::region_metaData_median_percent_methylation(coverage_all_chr, methylation_all_chr, region_metaData, col_oi = c(
    "Promoter", "exon", "intron",
    "CpG_island", "CpG_shelf", "CpG_shelf",
    "ChrSt_enhancer", "ChrSt_polycomb", "ChrSt_heterochromatin", "ChrSt_quies",
    "TE"
  ))

  print(patchwork::wrap_plots(plots, nrow = length(plots), ncol = 1))

  svglite(paste0("getcoverage_data_and_plots/RRBS_coverage_annotations_maxgap_",maxgap,".svg"), height = 12, fix_text_size = F)
  print(patchwork::wrap_plots(plots, nrow = length(plots), ncol = 1))
  dev.off()
}
