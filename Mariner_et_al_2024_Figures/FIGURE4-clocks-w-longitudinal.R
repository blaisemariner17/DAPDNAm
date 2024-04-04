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

pqlseq_res <- read_rds("../pqlseq_results/pqlseq_res_maxgap250.rds")
pqlseq_res <- pqlseq_res[pqlseq_res$converged == T,]
dim(pqlseq_res)

if (! "DAPDNAm" %in% installed.packages()){devtools::install_github("blaisemariner17/DAPDNAm", force = TRUE)}

# load in packages needed
library(glmnet)
library(tidyverse)

region_metaData <- read_rds("../metadata_regions/metaData_regions_maxgap_250.rds")

all(colnames(methylation_all_chr) == colnames(coverage_all_chr))
all(rownames(methylation_all_chr) == rownames(coverage_all_chr))

impute_mice_par <- function(perc_meth){
  for (row in rownames(perc_meth)){
    if(row == rownames(perc_meth)[1]){
      tempData <- mice::mice(t(perc_meth[c(1,2),]))
      perc_meth[c(paste(row)),] <- complete(tempData,1)[,paste(row)]
    }else{
      tempData <- mice::mice(t(perc_meth[c(paste(rownames(perc_meth)[1]),paste(row)),]))
      perc_meth[c(paste(row)),] <- complete(tempData,1)[,paste(row)]
    }
  }
  return(perc_meth)
}

read_in_perc_meth=T
if (read_in_perc_meth){
  perc_meth <- read_rds("../LongDat/perc_meth_imputed.rds")
} else {
  
  perc_meth <- methylation_all_chr / coverage_all_chr
  
  rm_row <- c()
  for (row_ in 1:nrow(perc_meth)) {
    row_of_interest <- perc_meth[row_,]
    if (sum(is.na(perc_meth[row_,])) == ncol(perc_meth)){rm_row <- append(rm_row, row_)}
  }
  perc_meth <- perc_meth[-c(rm_row),]

  chrs <- paste0("chr", c(1:38,"X"))
  perc_meth_list <- list()
  for(chr in chrs){
    chr_ <- paste0(chr,"_")
    perc_meth_list[[paste(chr)]] <- perc_meth[startsWith(rownames(perc_meth), chr_),]
  }
  perc_meth_list_imputed <- mclapply(
    perc_meth_list,
    impute_mice_par,
    mc.cores = 5
  )
  
  for (chr in chrs){
    if (chr == chrs[1]){
      perc_meth <- perc_meth_list_imputed[[paste(chr)]]
    }else{
      perc_meth <- rbind(perc_meth, perc_meth_list_imputed[[paste(chr)]])
    }
  }
  write_rds(perc_meth, "../LongDat/perc_meth_imputed.rds") 
  rm(perc_meth_list_imputed)
  rm(perc_meth_list)
}

# na_bysamples <- c()
# for (col in colnames(perc_meth)){
#   tot_ <- sum(is.na(perc_meth[,paste(col)]))
#   na_bysamples <- append(na_bysamples,tot_)
# }
# hist(na_bysamples)
# 
# na_byregion <- c()
# for (row_ in rownames(perc_meth)){
#   tot_ <- sum(is.na(perc_meth[paste(row_),]))
#   na_byregion <- append(na_byregion,tot_)
# }
# hist(na_byregion)

#### Promoter regions ####
perc_meth_prmtr <- perc_meth[rownames(perc_meth) %in% region_metaData$region[region_metaData$Promoter == 1],]
region_metaData_prmtr <- region_metaData[region_metaData$Promoter == 1,]
#### CpG Is regions ####
perc_meth_cpgisl <- perc_meth[rownames(perc_meth) %in% region_metaData$region[region_metaData$CpG_island == 1],]
region_metaData_cpgisl <- region_metaData[region_metaData$CpG_island == 1,]
#### TE regions #####
region_metaData_te <- region_metaData[region_metaData$TE == 1,]
perc_meth_te <- perc_meth[rownames(perc_meth) %in% region_metaData$region[region_metaData$TE == 1],]

set.seed(100)
alphs_oi <- list(0.05,0.1,0.15,0.2,0.25,
              0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,
              0.7,0.75,0.8,0.85,0.9,0.95,1)

time_start <- Sys.time()
print(time_start)
for (i in 1:10){
  print(i)
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
      # all_res_list[[paste(alph)]][["region_metaData"]] <- rbind(all_res_list[[paste(alph)]][["region_metaData"]],all_res_list_[[paste(alph)]][["region_metaData"]])
      all_res_list[[paste(alph)]][["clock_results"]] <- rbind(all_res_list[[paste(alph)]][["clock_results"]],all_res_list_[[paste(alph)]][["clock_results"]])
    }
  }
  print(Sys.time()- time_start)
}

for (alph in alphs_oi){
  print(alph)
  if (alph == alphs_oi[1]){
    all_res_clock <- all_res_list[[paste(alph)]][['clock_results']]
    all_res_metaData <- all_res_list[[paste(alph)]][['metaData']]
  } else {
    all_res_clock <- rbind(all_res_clock,all_res_list[[paste(alph)]][['clock_results']])
    all_res_metaData <- rbind(all_res_metaData,all_res_list[[paste(alph)]][['metaData']])
  }
}

# write_rds(x = all_res_metaData, file = "../LongDat/all_res_simulation_metaData.rds")
# write_rds(x = all_res_region_metaData, file = "../LongDat/all_res_simulation_region_metaData.rds")
write_rds(x = all_res_clock, file = "../LongDat/all_res_simulation_clock.rds")
write_rds(x = all_res_metaData, file = "../LongDat/all_res_simulation_metaData.rds")

all_res_clock <- read_rds("../LongDat/all_res_simulation_clock.rds")
all_res_metaData <- read_rds("../LongDat/all_res_simulation_metaData.rds")

#### plotting ####
all_res <- all_res_clock

all_res$label_gg <- gsub("_"," ",paste0(all_res$regions, " regions (", all_res$n_regions, ")"))

all_res_sum <- WLSplot::data_summary(all_res, groupnames = c("label_gg", "alpha"), varname = "MSE.")

colnames(all_res_sum) <- c("label_gg", "alpha", "MSE", "sd", "se")

plot1 <- ggplot(all_res_sum[all_res_sum$label_gg != "Random regions (12487)",], aes(x = alpha, y = MSE, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(ymin = MSE - se, ymax = MSE+se), width=0.02)+
  ggtitle(paste0("each dot represents ", i, " clock simulations"))+
  xlab("alpha")+
  ylab("Clock's mean standard error")+
  scale_color_manual(
    values = c("black","red", "darkgreen", "darkblue", "darkorange"),
  )+
  theme_blaise+
  labs(tag = "A")

plot1

for (region in unique(all_res_metaData$regions)){
  print(region)
  for (cohort in unique(all_res_metaData$Cohort)){
    print(cohort)
    duplicated_dogs <- duplicated(all_res_metaData$dog_id[all_res_metaData$Cohort == cohort &
                                                                                    all_res_metaData$regions == region])
    for (dog_id in unique(all_res_metaData$dog_id[duplicated_dogs])){
      mean_predicted <- mean(all_res_metaData$predicted[all_res_metaData$dog_id == dog_id & 
                                                          all_res_metaData$Cohort == cohort &
                                                          all_res_metaData$regions == region])
      all_res_metaData$predicted[all_res_metaData$dog_id == dog_id & 
                                   all_res_metaData$Cohort == cohort &
                                   all_res_metaData$regions == region] <- mean_predicted
    }
    all_res_metaData <- unique(all_res_metaData)
  }
}

plot_list <- list()
for (alph in c(0.25,0.5,0.75)){
  for (region in unique(all_res_metaData$regions)){
    if (region == "Random"){next}
    plot11 <- ggplot(all_res_metaData[all_res_metaData$regions == region & all_res_metaData$alpha == alph,],
                     aes(x = Age_at_sample, y = predicted, color = dog_id))+
      geom_abline(slope=1, color = "black")+
      geom_point()+
      geom_line()+
      ggtitle(paste0(gsub("_"," ",region), " regions \n alpha: ", alph))+
      ylab("Epigenetic age")+
      xlab("Age at sample collection")+
      theme_blaise + theme(legend.position = 'none')+# c(0.25, 0.82))+
      xlim(c(0,as.integer(max(all_res_metaData$Age_at_sample))+1))+
      ylim(c(0,as.integer(max(all_res_metaData$Age_at_sample))+1))+
      labs(tag = "")
    # print(plot11)
    plot_list[[paste0(region,alph)]] <- print(plot11)
  }
}

plot11 <- patchwork::wrap_plots(plot_list, nrow = 3, ncol = 4)
plot11

# if (list_res[['clock_results']]$MSE. <= min(list_res[['clock_results']]$MSE., all_res$MSE.[all_res$regions == "All"])){
#   plot1 <- ggplot(list_res[["metaData"]], aes(x = Age_at_sample, y = predicted, color = dog_id))+
#     geom_abline(slope=1, color = "black")+
#     geom_point()+
#     geom_line()+
#     ggtitle(paste0("All regions \n alpha: ", alph, " MSE: ", round(list_res[["clock_results"]]$MSE., digits = 3)))+
#     ylab("Epigenetic age")+
#     xlab("Age at sample collection")+
#     theme_blaise + theme(legend.position = 'none')+# c(0.25, 0.82))+
#     xlim(c(0,as.integer(max(list_res[["metaData"]]$Age_at_sample))+1))+
#     ylim(c(0,as.integer(max(list_res[["metaData"]]$Age_at_sample))+1))+
#     labs(tag = "A")
# }

# lets do a leave one out

for (breed_size in c("Giant","Small")){
  time_start <- Sys.time()
  print(time_start)
  all_res_list <- parallel::mclapply(alphs_oi, 
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
                                     unique(metaData$dog_id[metaData$Breed_size == breed_size]),
                                     mc.cores = 20
  )
  names(all_res_list) <- alphs_oi
  print(Sys.time()- time_start)
  
  for (alph in alphs_oi){
    print(alph)
    if (alph == alphs_oi[1]){
      all_res_metaData <- all_res_list[[paste(alph)]][['metaData']]
      all_res_clock <- all_res_list[[paste(alph)]][['clock_results']]
    } else {
      all_res_clock <- rbind(all_res_clock,all_res_list[[paste(alph)]][['clock_results']])
      all_res_metaData <- rbind(all_res_metaData,all_res_list[[paste(alph)]][['metaData']])
    }
  }
  all_res_clock$Breed_size <- breed_size
  if (breed_size == "Giant"){
    all_res_clock_final <- all_res_clock; all_res_metaData_final <- all_res_metaData
  } else {
    all_res_clock_final <- rbind(all_res_clock_final,all_res_clock); all_res_metaData_final <- rbind(all_res_metaData_final,all_res_metaData)
    }
}

#### plotting ####
write_rds(x = all_res_clock_final, file = "../LongDat/all_res_simulation_clock_breedsize.rds")
write_rds(x = all_res_metaData_final, file = "../LongDat/all_res_simulation_metadata_breedsize.rds")

all_res_clock_final <- read_rds("../LongDat/all_res_simulation_clock_breedsize.rds")
all_res_metaData_final <- read_rds("../LongDat/all_res_simulation_metadata_breedsize.rds")

all_res <- all_res_metaData_final

all_res$oi <- (all_res$predicted - all_res$Age_at_sample)
all_res$label_gg <- paste0(all_res$Breed_size, " dogs")

all_res_sum <- WLSplot::data_summary(all_res, groupnames = c("label_gg", "alpha", "regions"), varname = "oi")

plot2 <- ggplot(all_res_sum[all_res_sum$regions == "All",], aes(x = alpha, y = oi, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(ymin = oi - se, ymax = oi + se), width=0.02)+
  ggtitle(paste0(""))+
  xlab("alpha")+
  ylab("predicted age - actual age")+
  scale_color_manual(values = c("red", "black"))+
  theme_blaise+#theme(legend.position = 'none')+
  labs(tag = "B")

plot2

plot21 <- ggplot(all_res_sum[all_res_sum$regions == "TE",], aes(x = alpha, y = oi, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(ymin = oi - se, ymax = oi + se), width=0.02)+
  ggtitle(paste0("TE clock"))+
  xlab("alpha")+
  ylab("predicted age - actual age")+
  scale_color_manual(values = c("red", "black"))+
  theme_blaise+
  labs(tag = "C")
plot21

plot22 <- ggplot(all_res[all_res$regions == "All" & all_res$alpha == 0.5,], aes(x = label_gg, y = oi, color = label_gg, group = label_gg))+
  geom_quasirandom() +
  geom_boxplot(width = 0.2, color = "black", outliers = F)+
  # geom_errorbar(aes(ymin = oi - se, ymax = oi + se), width=0.02)+
  ggtitle(paste0(""))+
  xlab("")+
  ylab("predicted age - actual age")+
  scale_color_manual(values = c("red", "black"))+
  theme_blaise+  
  geom_signif(
    comparisons = list(c("Giant dogs", "Small dogs")),
    map_signif_level = T, textsize = 8, na.rm = T, test = 't.test', color = 'black', tip_length = 0 #step_increase = 0.05, 
  )+
  labs(tag = "C") 
plot22

t.test(all_res$oi[all_res$regions == "All" & all_res$alpha == 0.5 & all_res$Breed_size == "Small"], 
       all_res$oi[all_res$regions == "All" & all_res$alpha == 0.5 & all_res$Breed_size == "Giant"])


for (sex in c("Male", "Female")){
  time_start <- Sys.time()
  print(time_start)
  all_res_list <- parallel::mclapply(alphs_oi, 
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
                                     unique(metaData$dog_id[metaData$Sex == sex]),
                                     mc.cores = 20
  )
  names(all_res_list) <- alphs_oi
  print(Sys.time()- time_start)
  
  for (alph in alphs_oi){
    print(alph)
    if (alph == alphs_oi[1]){
      all_res_metaData <- all_res_list[[paste(alph)]][['metaData']]
      all_res_clock <- all_res_list[[paste(alph)]][['clock_results']]
    } else {
      all_res_clock <- rbind(all_res_clock,all_res_list[[paste(alph)]][['clock_results']])
      all_res_metaData <- rbind(all_res_metaData,all_res_list[[paste(alph)]][['metaData']])
    }
  }
  all_res_clock$Sex <- sex
  if (sex == "Male"){
    all_res_clock_final <- all_res_clock; all_res_metaData_final <- all_res_metaData
  } else {
    all_res_clock_final <- rbind(all_res_clock_final,all_res_clock); all_res_metaData_final <- rbind(all_res_metaData_final,all_res_metaData)
  }
}
write_rds(x = all_res_clock_final, file = "../LongDat/all_res_simulation_clock_sex.rds")
write_rds(x = all_res_metaData_final, file = "../LongDat/all_res_simulation_metadata_sex.rds")

all_res_clock_final <- read_rds("../LongDat/all_res_simulation_clock_sex.rds")
all_res_metaData_final <- read_rds("../LongDat/all_res_simulation_metadata_sex.rds")
#### plotting ####
all_res <- all_res_metaData_final

all_res$oi <- (all_res$predicted - all_res$Age_at_sample)
all_res$label_gg <- paste0(all_res$Sex, " dogs")

all_res_sum <- WLSplot::data_summary(all_res, groupnames = c("label_gg", "alpha", "regions"), varname = "oi")

plot4 <- ggplot(all_res_sum[all_res_sum$regions == "All",], aes(x = alpha, y = oi, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(ymin = oi - se, ymax = oi + se), width=0.02)+
  ggtitle(paste0(""))+
  xlab("alpha")+
  ylab("predicted age - actual age")+
  scale_color_manual(
    values = c("blue", "red")
    )+
  theme_blaise+#theme(legend.position = 'none')+
  labs(tag = "C")

plot4

plot41 <- ggplot(all_res_sum[all_res_sum$regions == "TE",], aes(x = alpha, y = oi, color = label_gg))+
  geom_point(size = 3)+
  geom_line()+
  geom_errorbar(aes(ymin = oi - se, ymax = oi + se), width=0.02)+
  ggtitle(paste0("TE clock"))+
  xlab("alpha")+
  ylab("predicted age - actual age")+
  scale_color_manual(
    values = c("blue", "red")
  )+
  theme_blaise+
  labs(tag = "E")

plot41

svglite("DAP_FIGURE4.svg", fix_text_size = F, height = 16, width = 12)
print(plot1 / plot2 / plot4)
dev.off()

svglite("DAP_FIGURE4_chrn_v_pred_allSims.svg", fix_text_size = F, height = 16, width = 16)
print(plot11)
dev.off()

save.image("FIG4.RData")
