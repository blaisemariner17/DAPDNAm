#### Blaise Mariner 
# get the differentially methylated regions in the demographics of intererst
## for questions contact bmarine2@asu.edu or blaisemariner17@gmail.com
## version of R 4.2.2 needed for bsseq as of the date below
## 2024-04-01
##
## abbreviations: oi = of interest; chr = chromosome; dap = dog aging project

library_list <- c(
  "ggplot2",
  "svglite",
  "ggExtra",
  "ggtext",
  "tidyverse",
  "corrplot",
  "glmnet",
  "PQLseq",
  "IRanges",
  "patchwork",
  "ggtext",
  "ggrepel",
  "ggbeeswarm",
  "ggsignif"
  
)
lapply(library_list, require, character.only = TRUE)
theme_blaise <- theme(plot.title.position = "plot", axis.text.x = element_text(angle=0),      plot.title = element_text(family = "sans", size = 14, hjust = 0.5, color="black", face='bold'),      plot.subtitle = element_text(family = "sans", size = 11, color="black"),      
                      axis.text = element_text(family = "sans", size = 18, color="black"),axis.title.y = element_markdown(family = "sans", size = 16),   
                      axis.title.x = element_markdown(family = "sans", size = 16),       panel.border = element_blank(),      axis.line = element_line(colour = "black", linewidth = 1),       axis.ticks = element_line(colour = "black", linewidth = 1),       legend.key.size = unit(1.5, 'cm'),      legend.key = element_rect(fill=NA),      legend.text = element_text(family = "sans", size = 20, hjust = 0),      legend.title = element_blank(),      legend.background = element_blank(),      legend.box.background = element_blank(),     panel.background = element_blank(),      panel.grid.major = element_line(colour = "black"),      panel.grid.minor = element_blank())+ removeGrid()

#1.5 hrs ish with 25 cores (1 managing) and 250 Gb mem (200 errors out) during peak hours
cores_ <- 24


# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

# metadata
p1 <- read_rds("../metadata_samples/P1-DAP-metaData-2400216.rds")
p2 <- read_rds("../metadata_samples/P2-DAP-metaData-240124.rds")
p3 <- read_rds("../metadata_samples/P3-DAP-metaData-240124.rds")

metaData <- rbind(p1,p2,p3)

maxgap = "250"

coverage_all_chr <- read_rds(paste0("../getcoverage_data_and_plots/longitudinal_maxGap_250_all_coverage_regions_oi.rds"))
methylation_all_chr <- read_rds(paste0("../getcoverage_data_and_plots/longitudinal_maxGap_250_all_methylation_regions_oi.rds"))
region_metaData <- readRDS(paste0("../metadata_regions/metaData_regions_maxgap_", maxgap,".rds"))

coverage_all_chr <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData$region,]
methylation_all_chr <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData$region,]

if (! "DAPDNAm" %in% installed.packages()){devtools::install_github("blaisemariner17/DAPDNAm", force = TRUE)}

# load in packages needed
library(glmnet)
library(tidyverse)

region_metaData <- read_rds("../metadata_regions/metaData_regions_maxgap_250.rds")

all(colnames(methylation_all_chr) == colnames(coverage_all_chr))
all(rownames(methylation_all_chr) == rownames(coverage_all_chr))

perc_meth <- methylation_all_chr / coverage_all_chr

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth)) {
  row_of_interest <- perc_meth[row_,]
  if (sum(is.na(perc_meth[row_,])) == ncol(perc_meth)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth[row_,])))
}
perc_meth <- perc_meth[-c(rm_row),]

#### Promoter regions ####
region_metaData_prmtr <- region_metaData[region_metaData$Promoter == 1,]
coverage_all_chr_prmtr <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData_prmtr$region,]
methylation_all_chr_prmtr <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData_prmtr$region,]
all(colnames(methylation_all_chr_prmtr) == colnames(coverage_all_chr_prmtr))
all(rownames(methylation_all_chr_prmtr) == rownames(coverage_all_chr_prmtr))
perc_meth_prmtr <- methylation_all_chr_prmtr / coverage_all_chr_prmtr

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth_prmtr)) {
  row_of_interest <- perc_meth_prmtr[row_,]
  if (sum(is.na(perc_meth_prmtr[row_,])) == ncol(perc_meth_prmtr)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth_prmtr[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth_prmtr[row_,])))
}
perc_meth_prmtr <- perc_meth_prmtr[-c(rm_row),]

#### CpG Is regions ####
region_metaData_cpgisl <- region_metaData[region_metaData$CpG_island == 1,]
coverage_all_chr_cpgisl <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData_cpgisl$region,]
methylation_all_chr_cpgisl <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData_cpgisl$region,]
all(colnames(methylation_all_chr_cpgisl) == colnames(coverage_all_chr_cpgisl))
all(rownames(methylation_all_chr_cpgisl) == rownames(coverage_all_chr_cpgisl))
perc_meth_cpgisl <- methylation_all_chr_cpgisl / coverage_all_chr_cpgisl

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth_cpgisl)) {
  row_of_interest <- perc_meth_cpgisl[row_,]
  if (sum(is.na(perc_meth_cpgisl[row_,])) == ncol(perc_meth_cpgisl)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth_cpgisl[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth_cpgisl[row_,])))
}
perc_meth_cpgisl <- perc_meth_cpgisl[-c(rm_row),]


#### TE regions #####
region_metaData_te <- region_metaData[region_metaData$TE == 1,]
coverage_all_chr_te <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData_te$region,]
methylation_all_chr_te <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData_te$region,]
all(colnames(methylation_all_chr_te) == colnames(coverage_all_chr_te))
all(rownames(methylation_all_chr_te) == rownames(coverage_all_chr_te))
perc_meth_te <- methylation_all_chr_te / coverage_all_chr_te

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth_te)) {
  row_of_interest <- perc_meth_te[row_,]
  if (sum(is.na(perc_meth_te[row_,])) == ncol(perc_meth_te)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth_te[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth_te[row_,])))
}
perc_meth_te <- perc_meth_te[-c(rm_row),]

set.seed(100)
alphs_oi <- list(0.05,0.1,0.15,0.2,0.25,
              0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,
              0.7,0.75,0.8,0.85,0.9,0.95,1)

time_start <- Sys.time()
print(time_start)
for (i in 1:25){
  all_res_list_ <- parallel::mclapply(alphs_oi, 
                                     DAPDNAm::clock_simulate,
                                     metaData, 
                                     region_metaData,
                                     perc_meth, 
                                     region_metaData_te,
                                     perc_meth_te,
                                     region_metaData_prmtr,
                                     perc_meth_prmtr,
                                     region_metaData_cpgisl,
                                     perc_meth_cpgisl,
                                     mc.cores = 20
  )
  names(all_res_list_) <- alphs_oi
  if (i == 1){
    all_res_list <- all_res_list_
  } else {
    for (alph in alphs_oi){
      all_res_list[[paste(alph)]][['metaData']] <- rbind(all_res_list[[paste(alph)]][['metaData']],all_res_list_[[paste(alph)]][['metaData']])
      all_res_list[[paste(alph)]][["region_metaData"]] <- rbind(all_res_list[[paste(alph)]][["region_metaData"]],all_res_list_[[paste(alph)]][["region_metaData"]])
      all_res_list[[paste(alph)]][["clock_results"]] <- rbind(all_res_list[[paste(alph)]][["clock_results"]],all_res_list_[[paste(alph)]][["clock_results"]])
    }
  }
}
print(time_start - Sys.time())

for (alph in alphs_oi){
  if (alph == alphs_oi[1]){
    all_res_metaData <- all_res_list[[paste(alph)]][['metaData']]
    all_res_region_metaData <- all_res_list[[paste(alph)]][['region_metaData']]
    all_res_clock <- all_res_list[[paste(alph)]][['clock_results']]
  } else {
    all_res_metaData <- rbind(all_res_metaData,all_res_list[[paste(alph)]][['metaData']])
    all_res_region_metaData <- rbind(all_res_region_metaData,all_res_list[[paste(alph)]][['region_metaData']])
    all_res_clock <- rbind(all_res_clock,all_res_list[[paste(alph)]][['clock_results']])
  }
}

write_rds(x = all_res_metaData, file = "../LongDat/all_res_simulation_metaData.rds")
write_rds(x = all_res_region_metaData, file = "../LongDat/all_res_simulation_region_metaData.rds")
write_rds(x = all_res_clock, file = "../LongDat/all_res_simulation_clock.rds")

#### plotting ####
all_res <- all_res_clock
all_res$MSE. <- all_res$MSE.

all_res$label_gg <- paste0(all_res$regions, " regions (", all_res$n_regions, ")")

all_res_sum <- WLSplot::data_summary(all_res, groupnames = c("label_gg", "alpha"), varname = "MSE.")

colnames(all_res_sum) <- c("label_gg", "alpha", "MSE", "sd", "se")

plot6 <- ggplot(all_res_sum, aes(x = alpha, y = MSE, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(ymin = MSE - se, ymax = MSE+se), width=0.02)+
  ggtitle(paste0("MSE calculations from ", i, " simulations"))+
  xlab("alpha")+
  ylab("Clock's mean standard error")+
  scale_color_manual(values = c("black", "darkred","darkorange", "darkblue", "darkgreen"))+
  theme_blaise+
  labs(tag = "D")
plot6

if (list_res[['clock_results']]$MSE. <= min(list_res[['clock_results']]$MSE., all_res$MSE.[all_res$regions == "All"])){
  plot1 <- ggplot(list_res[["metaData"]], aes(x = Age_at_sample, y = predicted, color = dog_id))+
    geom_abline(slope=1, color = "black")+
    geom_point()+
    geom_line()+
    ggtitle(paste0("All regions \n alpha: ", alph, " MSE: ", round(list_res[["clock_results"]]$MSE., digits = 3)))+
    ylab("Epigenetic age")+
    xlab("Age at sample collection")+
    theme_blaise + theme(legend.position = 'none')+# c(0.25, 0.82))+
    xlim(c(0,as.integer(max(list_res[["metaData"]]$Age_at_sample))+1))+
    ylim(c(0,as.integer(max(list_res[["metaData"]]$Age_at_sample))+1))+
    labs(tag = "A")
}
if (list_res_TE[['clock_results']]$MSE. <= min(list_res_TE[['clock_results']]$MSE., all_res$MSE.[all_res$regions == "TE"])){
  plot2 <- ggplot(list_res_TE[["metaData"]], aes(x = Age_at_sample, y = predicted, color = dog_id))+
    geom_abline(slope=1, color = "black")+
    geom_point()+
    geom_line()+
    ggtitle(paste0("TE regions \n alpha: ", alph, " MSE: ", round(list_res_TE[["clock_results"]]$MSE., digits = 3)))+
    ylab("Epigenetic age")+
    xlab("Age at sample collection")+
    theme_blaise + theme(legend.position = 'none')+#c(0.25, 0.82))+
    xlim(c(0,as.integer(max(list_res_TE[["metaData"]]$Age_at_sample))+1))+
    ylim(c(0,as.integer(max(list_res_TE[["metaData"]]$Age_at_sample))+1))+
    labs(tag = "B")
}
if (list_res_PR[['clock_results']]$MSE. <= min(list_res_PR[['clock_results']]$MSE., all_res$MSE.[all_res$regions == "Promoter"])){
  plot3 <- ggplot(list_res_PR[["metaData"]], aes(x = Age_at_sample, y = predicted, color = dog_id))+
    geom_abline(slope=1, color = "black")+
    geom_point()+
    geom_line()+
    ggtitle(paste0("Promoter regions \n alpha: ", alph, " MSE: ", round(list_res_PR[["clock_results"]]$MSE., digits = 3)))+
    ylab("Epigenetic age")+
    xlab("Age at sample collection")+
    theme_blaise + theme(legend.position = 'none')+#c(0.25, 0.82))+
    xlim(c(0,as.integer(max(list_res_PR[["metaData"]]$Age_at_sample))+1))+
    ylim(c(0,as.integer(max(list_res_PR[["metaData"]]$Age_at_sample))+1))+
    labs(tag = "C")
}
if (list_res_CpGIs[['clock_results']]$MSE. <= min(list_res_CpGIs[['clock_results']]$MSE., all_res$MSE.[all_res$regions == "CpG_island"])){
  plot4 <- ggplot(list_res_CpGIs[["metaData"]], aes(x = Age_at_sample, y = predicted, color = dog_id))+
    geom_abline(slope=1, color = "black")+
    geom_point()+
    geom_line()+
    ggtitle(paste0("CpG regions \n alpha: ", alph, " MSE: ", round(list_res_CpGIs[["clock_results"]]$MSE., digits = 3)))+
    ylab("Epigenetic age")+
    xlab("Age at sample collection")+
    theme_blaise + theme(legend.position = 'none')+#c(0.25, 0.82))+
    xlim(c(0,as.integer(max(list_res_CpGIs[["metaData"]]$Age_at_sample))+1))+
    ylim(c(0,as.integer(max(list_res_CpGIs[["metaData"]]$Age_at_sample))+1))+
    labs(tag = "C")
}
if (list_res_rand[['clock_results']]$MSE. <= min(list_res_rand[['clock_results']]$MSE., all_res$MSE.[all_res$regions == "Random"])){
  plot5 <- ggplot(list_res_rand[["metaData"]], aes(x = Age_at_sample, y = predicted, color = dog_id))+
    geom_abline(slope=1, color = "black")+
    geom_point()+
    geom_line()+
    ggtitle(paste0("Random regions \n alpha: ", alph, " MSE: ", round(list_res_rand[["clock_results"]]$MSE., digits = 3)))+
    ylab("Epigenetic age")+
    xlab("Age at sample collection")+
    theme_blaise + theme(legend.position = 'none')+#c(0.25, 0.82))+
    xlim(c(0,as.integer(max(list_res_rand[["metaData"]]$Age_at_sample))+1))+
    ylim(c(0,as.integer(max(list_res_rand[["metaData"]]$Age_at_sample))+1))+
    labs(tag = "C")
}


svglite("DAP_FIGURE4.svg", fix_text_size=F, width = 14, height = 14)
(plot1 + plot2 +plot3 ) / plot6
dev.off()

#what if I train with precision_1 samples and test with 2,3?#

for (alph in alphs_oi){
  print(alph)
  list_res <- DAPDNAm::build_clock(alph,
                                   metaData,
                                   region_metaData,
                                   perc_meth,
                                   test_samples = colnames(perc_meth) %in% metaData$lid_pid[metaData$Cohort != "precision_1"]
  )
  
  list_res_TE <- DAPDNAm::build_clock(alph,
                                      metaData,
                                      region_metaData_te,
                                      perc_meth_te,
                                      test_samples = colnames(perc_meth) %in% metaData$lid_pid[metaData$Cohort != "precision_1"]
  )
  
  list_res_PR <- DAPDNAm::build_clock(alph,
                                      metaData,
                                      region_metaData_prmtr,
                                      perc_meth_prmtr,
                                      test_samples = colnames(perc_meth) %in% metaData$lid_pid[metaData$Cohort != "precision_1"]
  )
  
  list_res_CpGIs <- DAPDNAm::build_clock(alph,
                                         metaData,
                                         region_metaData_cpgisl,
                                         perc_meth_cpgisl,
                                         test_samples = colnames(perc_meth) %in% metaData$lid_pid[metaData$Cohort != "precision_1"]
  )
  
  list_res_rand <- DAPDNAm::build_clock(alph,
                                        metaData,
                                        region_metaData_rand,
                                        perc_meth_rand,
                                        test_samples = colnames(perc_meth) %in% metaData$lid_pid[metaData$Cohort != "precision_1"]
  )
  
  list_res[['clock_results']]$regions <- "All"
  list_res_TE[['clock_results']]$regions <- "TE"
  list_res_PR[['clock_results']]$regions <- "Promoter"
  list_res_CpGIs[['clock_results']]$regions <- "CpG_island"
  list_res_rand[['clock_results']]$regions <- "Random"
  
  if (alph == alphs_oi[1]){
    all_res <- rbind(list_res[['clock_results']], list_res_TE[['clock_results']],list_res_PR[['clock_results']], list_res_CpGIs[['clock_results']],list_res_rand[['clock_results']])
  }else{
    all_res <- rbind(all_res, list_res[['clock_results']], list_res_TE[['clock_results']],list_res_PR[['clock_results']], list_res_CpGIs[['clock_results']],list_res_rand[['clock_results']])
  }
  
  
  save_res <- all_res
}

#### plotting ####
all_res <- save_res
all_res <- all_res[order(all_res$regions),]

colnames(all_res) <- c("alpha", "MSE", "n_regions", "training_samples", "testing_samples", "regions")
all_res <- all_res[,c("regions", "n_regions", "alpha", "MSE", "training_samples", "testing_samples")]

all_res$label_gg <- paste0(all_res$regions, " regions (", all_res$n_regions, ")")

all_res<-all_res[all_res$alpha>0,]

plot6 <- ggplot(all_res, aes(x = alpha, y = MSE, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  ggtitle("p1 training")+
  xlab("alpha")+
  ylab("Clock's mean standard error")+
  scale_color_manual(values = c("black", "darkred","darkorange", "darkblue", "darkgreen"))+
  theme_blaise+
  labs(tag = "D")
plot6


# what if I only look at the dogs below age 10? #

# metadata
p1 <- read_rds("../metadata_samples/P1-DAP-metaData-2400216.rds")
p2 <- read_rds("../metadata_samples/P2-DAP-metaData-240124.rds")
p3 <- read_rds("../metadata_samples/P3-DAP-metaData-240124.rds")

metaData <- rbind(p1,p2,p3)
metaData <- metaData[metaData$Age_at_sample < 10,]

maxgap = "250"

coverage_all_chr <- read_rds(paste0("../getcoverage_data_and_plots/longitudinal_maxGap_250_all_coverage_regions_oi.rds"))
methylation_all_chr <- read_rds(paste0("../getcoverage_data_and_plots/longitudinal_maxGap_250_all_methylation_regions_oi.rds"))
region_metaData <- readRDS(paste0("../metadata_regions/metaData_regions_maxgap_", maxgap,".rds"))

coverage_all_chr <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData$region,]
methylation_all_chr <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData$region,]

if (! "DAPDNAm" %in% installed.packages()){devtools::install_github("blaisemariner17/DAPDNAm", force = TRUE)}

# load in packages needed
library(glmnet)
library(tidyverse)

region_metaData <- read_rds("../metadata_regions/metaData_regions_maxgap_250.rds")

all(colnames(methylation_all_chr) == colnames(coverage_all_chr))
all(rownames(methylation_all_chr) == rownames(coverage_all_chr))

perc_meth <- methylation_all_chr / coverage_all_chr

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth)) {
  row_of_interest <- perc_meth[row_,]
  if (sum(is.na(perc_meth[row_,])) == ncol(perc_meth)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth[row_,])))
}
perc_meth <- perc_meth[-c(rm_row),]

#### Promoter regions ####
region_metaData_prmtr <- region_metaData[region_metaData$Promoter == 1,]
coverage_all_chr_prmtr <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData_prmtr$region,]
methylation_all_chr_prmtr <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData_prmtr$region,]
all(colnames(methylation_all_chr_prmtr) == colnames(coverage_all_chr_prmtr))
all(rownames(methylation_all_chr_prmtr) == rownames(coverage_all_chr_prmtr))
perc_meth_prmtr <- methylation_all_chr_prmtr / coverage_all_chr_prmtr

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth_prmtr)) {
  row_of_interest <- perc_meth_prmtr[row_,]
  if (sum(is.na(perc_meth_prmtr[row_,])) == ncol(perc_meth_prmtr)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth_prmtr[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth_prmtr[row_,])))
}
perc_meth_prmtr <- perc_meth_prmtr[-c(rm_row),]

#### CpG Is regions ####
region_metaData_cpgisl <- region_metaData[region_metaData$CpG_island == 1,]
coverage_all_chr_cpgisl <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData_cpgisl$region,]
methylation_all_chr_cpgisl <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData_cpgisl$region,]
all(colnames(methylation_all_chr_cpgisl) == colnames(coverage_all_chr_cpgisl))
all(rownames(methylation_all_chr_cpgisl) == rownames(coverage_all_chr_cpgisl))
perc_meth_cpgisl <- methylation_all_chr_cpgisl / coverage_all_chr_cpgisl

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth_cpgisl)) {
  row_of_interest <- perc_meth_cpgisl[row_,]
  if (sum(is.na(perc_meth_cpgisl[row_,])) == ncol(perc_meth_cpgisl)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth_cpgisl[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth_cpgisl[row_,])))
}
perc_meth_cpgisl <- perc_meth_cpgisl[-c(rm_row),]

#### TE regions #####
region_metaData_te <- region_metaData[region_metaData$TE == 1,]
coverage_all_chr_te <- coverage_all_chr[rownames(coverage_all_chr) %in% region_metaData_te$region,]
methylation_all_chr_te <- methylation_all_chr[rownames(methylation_all_chr) %in% region_metaData_te$region,]
all(colnames(methylation_all_chr_te) == colnames(coverage_all_chr_te))
all(rownames(methylation_all_chr_te) == rownames(coverage_all_chr_te))
perc_meth_te <- methylation_all_chr_te / coverage_all_chr_te

# tempData <- mice::mice(t(epi))
# completedData <- mice::complete(tempData,1)
rm_row <- c()
for (row_ in 1:nrow(perc_meth_te)) {
  row_of_interest <- perc_meth_te[row_,]
  if (sum(is.na(perc_meth_te[row_,])) == ncol(perc_meth_te)){rm_row <- append(rm_row, row_)}
  mean_row <- mean(row_of_interest, na.rm = T)
  row_of_interest[is.na(row_of_interest)] <- mean_row
  perc_meth_te[row_,] <- row_of_interest
  # print(sum(is.na(perc_meth_te[row_,])))
}
perc_meth_te <- perc_meth_te[-c(rm_row),]

#### random regions #####
set.seed(100)
sampling <- sample(rownames(perc_meth), as.integer(mean(nrow(perc_meth_te), nrow(perc_meth_prmtr), nrow(perc_meth_cpgisl))))
perc_meth_rand <- perc_meth[rownames(perc_meth) %in% sampling,]
region_metaData_rand <- region_metaData[region_metaData$region %in% sampling,]

sampling <- sample(unique(metaData$dog_id), 90)

for (alph in alphs_oi){
  print(alph)
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
  
  list_res[['clock_results']]$regions <- "All"
  list_res_TE[['clock_results']]$regions <- "TE"
  list_res_PR[['clock_results']]$regions <- "Promoter"
  list_res_CpGIs[['clock_results']]$regions <- "CpG_island"
  list_res_rand[['clock_results']]$regions <- "Random"
  
  if (alph == alphs_oi[1]){
    all_res <- rbind(list_res[['clock_results']], list_res_TE[['clock_results']],list_res_PR[['clock_results']], list_res_CpGIs[['clock_results']],list_res_rand[['clock_results']])
  }else{
    all_res <- rbind(all_res, list_res[['clock_results']], list_res_TE[['clock_results']],list_res_PR[['clock_results']], list_res_CpGIs[['clock_results']],list_res_rand[['clock_results']])
  }
  
  save_res <- all_res
}

#### plotting ####
all_res <- save_res
all_res <- all_res[order(all_res$regions),]

colnames(all_res) <- c("alpha", "MSE", "n_regions", "training_samples", "testing_samples", "regions")
all_res <- all_res[,c("regions", "n_regions", "alpha", "MSE", "training_samples", "testing_samples")]

all_res$label_gg <- paste0(all_res$regions, " regions (", all_res$n_regions, ")")

all_res<-all_res[all_res$alpha>0,]

plot6 <- ggplot(all_res, aes(x = alpha, y = MSE, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  ggtitle("just dogs samples w age < 10")+
  xlab("alpha")+
  ylab("Clock's mean standard error")+
  scale_color_manual(values = c("black", "darkred","darkorange", "darkblue", "darkgreen"))+
  theme_blaise+
  labs(tag = "D")
plot6