#### Blaise Mariner 
#### this is the script used to generate the sample metadata for all RRBS dogs
## 2024-05-09

library_list <- c(
  "Biostrings",
  "parallel",
  "stringr"
)

lapply(library_list, require, character.only = TRUE)

# this sets the working directory to this script's path
if (isRStudio <- Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
print(getwd())

if (dir.exists("metadata_samples") == FALSE) {dir.create("metadata_samples")}

###########################################################################################################################
###########################################################################################################################
############################## read in all the  necessary information #####################################################
###########################################################################################################################
###########################################################################################################################

# # save out the lids if you need to pull from shiny
# dap <- readRDS("/scratch/nsnyderm/dap_rrbs/bismarkBSseq.rds")
# sampleNames(dap) <- paste(sub(".*.(lid_\\d+)_.*", "\\1", sampleNames(dap)), 
#                           sub(".*_(pid_\\d+).*", "\\1", sampleNames(dap)), 
#                           sep = "_")
# lids <- gsub("(lid_\\d+)_.*", "\\1", sampleNames(dap))
# write.csv(as.data.frame(lids), file = 'metadata_samples/lids.csv')

# statsQC <- read.csv("../colData_metaData_final_240216/dap_rrbsQC_stats-240216.csv")
statsQC <- read.csv("../colData_metaData_final_240216/dap_rrbsQC_stats.txt", sep = " ")

#### Load in the DAP, Shiny, and Clinical metadata
#DAP data
load("../colData_metaData_final_240216/DAP_2023_DogOverview_v1.0.RData")
for (column_ in 1:ncol(DogOverview)) {
  DogOverview[,column_] <- as.character(haven::as_factor(DogOverview[,column_]))
}

#shiny data, manually scraped at https://shiny.geladaresearch.org/
shiny_metaData <- read.csv("metadata_samples/shiny_240508.csv")
shiny_metaData$sid <- as.character(shiny_metaData$dap_sample_id)
shiny_metaData$dap_elab_cohort <- gsub(" ", "_",shiny_metaData$dap_elab_cohort)
shiny_metaData$dap_elab_cohort <- tolower(shiny_metaData$dap_elab_cohort)

#clinical data
clinical_metaData <- read.csv("../colData_metaData_final_240216/DAP_2023_SamplesResults_Metadata_v1.0.csv")
clinical_metaData$dog_id <- as.character(clinical_metaData$dog_id)
clinical_metaData <- clinical_metaData[clinical_metaData$Sample_Type %in% c("RRBS"),]

#most up-to-date precision dogs
pauls_precision <- read.csv("../colData_metaData_final_240216/precision240301-fromPaul.csv")
pauls_precision$dog_id <- as.character(pauls_precision$dog_id)

#### FILTER OUT THE LOW READS ####
lowreads <- statsQC[statsQC$reads < 1e6,]
rownames(lowreads) <- lowreads$lid_pid
write.csv("metadata_samples/LOWREADS_FILTERED_OUT.csv", x = lowreads)

statsQC <- statsQC[statsQC$reads > 1e6,]
statsQC$lid <- gsub("_PID_.*", "", statsQC$lid_pid)
statsQC$pid <- gsub(".*_PI", "PI", statsQC$lid_pid)

###########################################################################################################################
###########################################################################################################################
#############  First step lets make a data frame with the lid_pid, lid, pid, and reads ####################################
###########################################################################################################################
###########################################################################################################################

# remove metaData if it exists
if(exists("metaData") == T){rm(metaData)}
statsQC$Perc_unique <- statsQC$unique / statsQC$reads
statsQC$perc_meth <- statsQC$meth_cpg / (statsQC$meth_cpg + statsQC$unmeth_cpg)
metaData <- statsQC[,c("lid_pid", "lid", "pid", "reads", "unique", "Perc_unique", "chrX_ratio", "perc_meth")]
rownames(metaData) <- metaData$lid_pid

#add dap sample id, since some lid's are run twice we have to do a for loop here
for (lid in unique(metaData$lid)) {
  metaData$sid[metaData$lid == lid] <- shiny_metaData$sid[shiny_metaData$lid == lid]
  metaData$rdid[metaData$lid == lid] <- shiny_metaData$parent_id[shiny_metaData$lid == lid]
  metaData$prep_date[metaData$lid == lid] <- shiny_metaData$prep_date[shiny_metaData$lid == lid]
  metaData$shiny_dog_id[metaData$lid == lid] <- as.character(shiny_metaData$dap_dog_id[shiny_metaData$lid == lid])
}

###########################################################################################################################
###########################################################################################################################
############################## some elab formatting is fucked up so this is needed ########################################
###########################################################################################################################
###########################################################################################################################

# #note any NA values for dap sample ids
# # sum(is.na(metaData$sid))
# write.csv(x = metaData[is.na(metaData$sid),], file = "metadata_samples/na_shiny_metadata_searching/na_sample_id_shiny.csv")
# # Then you're going to use the eLab metadata shiny app to get the RDIDs and save that file and upload it here:
# elab_wrong_format_rdid_file <- read.csv("metadata_samples/na_shiny_metadata_searching/rdid_na_samples.csv")
# 
# #then we're going to scrape through with the appropriate id (this is the integer associated with DAP-####-#####D0p number, idk how else to get it besides this way)
# library(httr)
# library(jsonlite)
# # load R script that loads token
# auth_token="9da7eacfb173556e6d12fdd6b469a742"
# auth_token2="49644bccf71d81754be8f37a164b3125"
# get_sid <- function(id) {
#   sampleID=gsub(" ","%20",id)
#   get_meta = paste0("https://us.elabjournal.com/api/v1/samples/",sampleID,"/meta")
#   call_api <- httr::GET(url = get_meta, 
#                         httr::add_headers(.headers = c('Authorization'= auth_token2, 'Accept' = 'text/json')))
#   out = parse_json(httr::content(call_api))
#   return(out$data[[1]]$value)
# }
# elab_wrong_format_rdid_file$Source <- gsub(".*\\|", "", elab_wrong_format_rdid_file$Source)
# for (lid in metaData$lid){
#   rdid <- unique(metaData$rdid[metaData$lid == lid])
#   if (rdid %in% elab_wrong_format_rdid_file$ID){
#     DAP_idk_wtf_id_this_is <- elab_wrong_format_rdid_file$Source[elab_wrong_format_rdid_file$ID == rdid]
#     correct_dap_sid <- get_sid(DAP_idk_wtf_id_this_is)
#     print(correct_dap_sid)
#     metaData$sid[metaData$lid == lid] <- correct_dap_sid 
#   }
# }
# 
# #this should now equal zero :)
# sum(is.na(metaData$sid))

###########################################################################################################################
###########################################################################################################################
########################## now we can add the dog_ids from the clinical ###################################################
###########################################################################################################################
###########################################################################################################################

#most of the dog_ids are going to be in the clinical
for (sid in metaData$sid) {
  clinical_dog_id <- unique(clinical_metaData$dog_id[clinical_metaData$DAP_Sample_ID == sid])
  if (length(clinical_dog_id) == 0) {clinical_dog_id = NA}
  metaData$clinical_dog_id[metaData$sid == sid] <- clinical_dog_id
  
  shiny_dog_id <- unique(metaData$shiny_dog_id[metaData$sid == sid])
  if (length(shiny_dog_id) == 0){shiny_dog_id = NA}
  
  if (!is.na(shiny_dog_id) & !is.na(clinical_dog_id)){
    if (shiny_dog_id == clinical_dog_id){
      metaData$dog_id[metaData$sid == sid] <- shiny_dog_id
    } else {
      metaData$dog_id[metaData$sid == sid] <- "contradiction"
      print("uh oh")
    }
  }
  if (is.na(shiny_dog_id) & !is.na(clinical_dog_id)){
    metaData$dog_id[metaData$sid == sid] <- clinical_dog_id
  }
  
  if (!is.na(shiny_dog_id) & is.na(clinical_dog_id)){
    metaData$dog_id[metaData$sid == sid] <- shiny_dog_id
  }
}
write.csv(file = "metadata_samples/na_dog_id_clinical_and_eLab.csv", metaData[is.na(metaData$shiny_dog_id) & is.na(metaData$clinical_dog_id),])

###########################################################################################################################
###########################################################################################################################
########################## now we can add the cohort information####### ###################################################
###########################################################################################################################
###########################################################################################################################

#triad cohort information is going to be from the shiny_metaData
for (sid in metaData$sid){
  cohort <- unique(na.omit(shiny_metaData$dap_elab_cohort[shiny_metaData$sid == sid]))
  if(length(cohort) == 0){next}
  if (grepl("triad", cohort)) {
    metaData$Cohort_group[metaData$sid == sid] <- "Triad"
    metaData$Cohort[metaData$sid == sid] <- cohort
  }
}
# precision information is from paul's metadata
for (sid in metaData$sid){
  dog_id <- unique(na.omit(metaData$dog_id[metaData$sid ==sid]))
  if (length(dog_id) == 0){next}
  if (dog_id %in% pauls_precision$dog_id){
    metaData$Cohort_group[metaData$sid == sid] <- "Precision"
    cohort <- clinical_metaData$Sample_Year[clinical_metaData$DAP_Sample_ID == sid]
    if (length(cohort) == 0){cohort = "need updated clincial_metaData file"}
    metaData$Cohort[metaData$sid == sid] <- cohort
  }
}

###########################################################################################################################
###########################################################################################################################
##################### now we can add demographic information from DO and clinical##########################################
###########################################################################################################################
###########################################################################################################################

not_in_DogOv <- c()
not_in_clinical <- c()
for (sid in unique(metaData$sid)){
  dog_id <- unique(metaData$dog_id[metaData$sid == sid])
  if (length(dog_id) == 0) {next}
  if (is.na(dog_id)) {next}
  if(!dog_id %in% DogOverview$dog_id){not_in_DogOv <- append(not_in_DogOv, dog_id);next}
  
  #DO dems
  LifeStage_Class_at_HLES <- unique(na.omit(DogOverview$LifeStage_Class_at_HLES[DogOverview$dog_id == dog_id]));metaData$LifeStage_Class_at_HLES[metaData$sid == sid] <- LifeStage_Class_at_HLES
  Estimated_DOB <- DogOverview$Estimated_DOB[DogOverview$dog_id == dog_id];metaData$Estimated_DOB[metaData$sid == sid] <- Estimated_DOB
  Breed_Status <- DogOverview$Breed_Status[DogOverview$dog_id == dog_id];metaData$Breed_Status[metaData$sid == sid] <- Breed_Status
  breed <- DogOverview$Breed[DogOverview$dog_id == dog_id]; metaData$Breed[metaData$sid == sid] <- breed
  Sex_class <- DogOverview$Sex_class_at_HLES[DogOverview$dog_id == dog_id]; metaData$Sex_class[metaData$sid == sid] <- Sex_class
  Breed_Size_Class_at_HLES <- DogOverview$Breed_Size_Class_at_HLES[DogOverview$dog_id == dog_id]; metaData$Breed_Size_Class_at_HLES[metaData$sid == sid] <- Breed_Size_Class_at_HLES
  
  #clinical info
  if (!sid %in% clinical_metaData$DAP_Sample_ID){not_in_clinical <- append(not_in_clinical, sid);next}
  Sample_taken <- unique(na.omit(clinical_metaData$Sample_Collection_DateTime[clinical_metaData$DAP_Sample_ID == sid]))
  if (length(Sample_taken) > 0) {if (grepl("/",Sample_taken)){split <- str_split(Sample_taken, "/");month <- split[[1]][1];day <- split[[1]][2];split <- str_split(split[[1]][3], " ");year <- split[[1]][1];time <- split[[1]][2];Sample_taken <- paste0(year, "-", month, "-", day, " ", time)}}
  metaData$Sample_taken[metaData$sid == sid] <- Sample_taken
  
  #the weight dem is a pain in my ass but we got it now
  weight_at_sample_collection <- clinical_metaData$Sample_Dog_Weight[clinical_metaData$DAP_Sample_ID == sid]
  weight_units <- (clinical_metaData$Sample_Dog_Weight_Units[clinical_metaData$DAP_Sample_ID == sid])
  if(length(weight_units) > 0){if (is.na(weight_units)==F){if (weight_units == "lb") {weight_at_sample_collection = weight_at_sample_collection * 0.45359237;weight_units = "kg"}}}
  if (length(weight_at_sample_collection) == 0) {weight_at_sample_collection <- NA}
  if (is.na(weight_at_sample_collection)){
    weight_at_sample_collection <- DogOverview$Weight_Class_5KGBin_at_HLES[DogOverview$dog_id == dog_id]
    weight_at_sample_collection <- substr(weight_at_sample_collection, start = 1, stop = 2)
    weight_at_sample_collection <- gsub("-","",weight_at_sample_collection)
    weight_at_sample_collection <- as.numeric(weight_at_sample_collection) + 5
    weight_units = "kg"
  }
  metaData$weight[metaData$sid == sid] <- weight_at_sample_collection
  metaData$weight_units[metaData$sid == sid] <- weight_units
}

# not_in_DogOv #these are likely triad dogs
write.csv(not_in_DogOv, file = "metadata_samples/DAP_dog_ids_notin_DogOverview.csv")
#these are likely new samples as they are in the dog overview
oi <- metaData[metaData$sid %in% not_in_clinical,]
write.csv(oi, file = "metadata_samples/DAP_SIDs_notin_clinical_metadata.csv")

metaData$Sex <- "unknown"
metaData$Sex[grepl(pattern = "Male", metaData$Sex_class)] <- "Male"
metaData$Sex[grepl(pattern = "Female", metaData$Sex_class)] <- "Female"

metaData$Breed[grepl(pattern = "nknown", x = metaData$Breed)] <- "unknown"

metaData$fixed <- "unknown"
metaData$fixed[grepl(pattern = "spayed", x = metaData$Sex_class)] <- "Spayed"
metaData$fixed[grepl(pattern = "neutered", x = metaData$Sex_class)] <- "Neutered"
metaData$fixed[grepl(pattern = "intact", x = metaData$Sex_class)] <- "Intact"

metaData$Fixed_bool[metaData$fixed == "Intact"] <- 0
metaData$Fixed_bool[metaData$fixed == "Spayed"] <- 1
metaData$Fixed_bool[metaData$fixed == "Neutered"] <- 1

metaData$Sex_bool[metaData$Sex == "Female"] <- 1
metaData$Sex_bool[metaData$Sex == "Male"] <- 0

metaData$Breed_Status_bool[metaData$Breed_Status == "Purebred"] <- 0
metaData$Breed_Status_bool[metaData$Breed_Status == "Mixed breed"] <- 1

metaData$Breed_size[grepl("Toy",metaData$Breed_Size_Class_at_HLES)] <- "Small"
metaData$Breed_size[grepl("Medium",metaData$Breed_Size_Class_at_HLES)] <- "Medium"
metaData$Breed_size[grepl("Standard",metaData$Breed_Size_Class_at_HLES)] <- "Standard"
metaData$Breed_size[grepl("Large",metaData$Breed_Size_Class_at_HLES)] <- "Large"
metaData$Breed_size[grepl("Giant",metaData$Breed_Size_Class_at_HLES)] <- "Giant"

################################################################################# DEMOGRAHIC OR CLINICAL CORRECTIONS: 
metaData$weight[metaData$dog_id == "94533"] <- 78 / 2.2
metaData$weight[metaData$dog_id == "118580"] <- 46.4 / 2.2
metaData$weight_units[metaData$dog_id == "94533"] <- "kg"
metaData$weight_units[metaData$dog_id == "94533"] <- "kg"
metaData$Breed[metaData$dog_id == "118580"] <- "Boxer"
metaData$LifeStage_Class_at_HLES[metaData$dog_id %in% c("118580", "94533")] <- "Puppy"
metaData$Breed_Size_Class_at_HLES[metaData$dog_id %in% c("94533")] <- "Large non-AKC and mixed breed dogs"
metaData$Breed_Size_Class_at_HLES[metaData$dog_id %in% c("118580")] <- "Standard non-AKC breed dogs (20 - 29.9 kg)"
metaData$Estimated_DOB[metaData$dog_id == "60812"] <- "2005-11-04"
metaData$Estimated_DOB[metaData$dog_id == "63753"] <- "2013-01-01"
metaData$Estimated_DOB[metaData$dog_id == "83095"] <- "2016-01-10"
metaData$Estimated_DOB[metaData$dog_id == "99570"] <- "2017-01-31"
#practice dogs 78968 & 79574 & 300036
metaData <- metaData[!metaData$dog_id %in% c('78968','79574', '300036'),]

###########################################################################################################################
###########################################################################################################################
##################### Age at sample collection calculation ################################################################
###########################################################################################################################
###########################################################################################################################

metaData$Estimated_DOB[is.na(metaData$Estimated_DOB)] <- "unknown"
metaData$Age_at_sample <- "unknown"
for (sid in metaData$sid[metaData$Estimated_DOB != "unknown"]){
  dob <- paste(metaData$Estimated_DOB, "00:00:00")
  metaData$Age_at_sample[metaData$sid == sid] <- as.numeric(as.POSIXct(metaData$Sample_taken[metaData$sid == sid]) - as.POSIXct(metaData$Estimated_DOB[metaData$sid == sid])) / 365.25
}

###########################################################################################################################
###########################################################################################################################
#####################      Predicted height    ###########################################################################
###########################################################################################################################
###########################################################################################################################

#predicted height
predicted_height <- read.csv("../colData_metaData_final_240216/precision_dogs_N-837_bodysize.phenotypes-240214.csv")
predicted_height$dog_id <- as.character(predicted_height$id)

for (sid in metaData$sid) {
  dog_id <- unique(na.omit(metaData$dog_id[metaData$sid == sid]))
  if (length(dog_id) == 0){next}
  height <- predicted_height$prediction[predicted_height$dog_id == dog_id]
  if (length(height) != 0){
    weight <- as.numeric(metaData$weight[metaData$sid == sid]); weight_over_height <- weight / height
    metaData$weight_over_height[metaData$sid == sid] <- weight_over_height
  }
}

###########################################################################################################################
###########################################################################################################################
#####################           SAMPLE SWAPS    ###########################################################################
###########################################################################################################################
###########################################################################################################################

genotyping_file <- readRDS("metadata_samples/dap_check_lid_res.rds")
genotyping_file$lid_pid <- gsub(".match", "", genotyping_file$lid_pid)
genotyping_file$pid <- gsub(".*_P", "P", genotyping_file$lid_pid)
genotyping_file <- genotyping_file[genotyping_file$lid_pid %in% metaData$lid_pid,]
genotyping_file <- genotyping_file[order(genotyping_file$lid_pid),]
metaData <- metaData[order(metaData$lid_pid),]
all(metaData$lid_pid == genotyping_file$lid_pid)

genotyping_res <- cbind(metaData[,c("lid_pid", "lid", "pid", "rdid", "sid", "dog_id")], genotyping_file[,c(4,5,6,1,2,3)])
genotyping_res$tophit <- as.character(genotyping_res$tophit)
genotyping_res$secondhit <- as.character(genotyping_res$secondhit)
genotyping_res$thirdhit <- as.character(genotyping_res$thirdhit)
mismatches <- genotyping_res
for (lid_pid in unique(mismatches$lid_pid)){
  if (is.na(unique(mismatches$dog_id[mismatches$lid_pid == lid_pid]))) {next}
  if (is.na(unique(mismatches$tophit[mismatches$lid_pid == lid_pid]))) {next}

  if (unique(mismatches$dog_id[mismatches$lid_pid == lid_pid]) == unique(mismatches$tophit[mismatches$lid_pid == lid_pid])){
    mismatches <- mismatches[mismatches$lid_pid != lid_pid,]
  }
}
mismatches <- mismatches[!is.na(mismatches$lid),]
mismatches <- mismatches[order(mismatches$pid),]
metaData$no_lid <- as.numeric(gsub("LID_", "", metaData$lid))
metaData$no_pid <- as.numeric(gsub("PID_", "", metaData$pid))
pid_oi <- metaData[metaData$pid == "PID_10558" , ]
pid_oi <- pid_oi[order(pid_oi$no_lid),]
pid_oi[,c("pid","lid","dog_id")]
#remove the identified swaps so far
mismatches <- mismatches[!rownames(mismatches) %in% c("LID_1070233_PID_10558", "LID_1070236_PID_10558", 
                                                      "LID_1071350_PID_10558", "LID_1070234_PID_10558"
                                                      ),]

pool_oi <- mismatches[mismatches$pid == "PID_10558"]
dog_id_oi = "97916"
pool_oi[pool_oi$dog_id == dog_id_oi | pool_oi$tophit == dog_id_oi| pool_oi$secondhit == dog_id_oi | pool_oi$thirdhit ==dog_id_oi,]
#swaps identified 5/10/24 blm
metaData$lid_pid[metaData$lid_pid == "LID_1070233_PID_10558" & metaData$lid == "LID_1070233"] <- "_LID_1070236_PID_10558"
metaData$lid_pid[metaData$lid_pid == "LID_1070236_PID_10558" & metaData$lid == "LID_1070236"] <- "_LID_1070233_PID_10558"
metaData$lid_pid[metaData$lid_pid == "LID_1070234_PID_10558" & metaData$lid == "LID_1070234"] <- "_LID_1071350_PID_10558"
metaData$lid_pid[metaData$lid_pid == "LID_1071350_PID_10558" & metaData$lid == "LID_1071350"] <- "_LID_1071334_PID_10558"
# metaData$lid_pid[metaData$lid_pid == "LID_1070232_PID_10558" & metaData$lid == "LID_1070232"] <- "_LID_1070237_PID_10558"
# metaData$lid_pid[metaData$lid_pid == "LID_1070230_PID_10558" & metaData$lid == "LID_1070230"] <- "_LID_1070223_PID_10558"

#swaps in pid 10559
pool_oi <- mismatches[mismatches$pid == "PID_10559"
           ,]
for (dog_id in pool_oi$dog_id){
  mismatch_hit <- pool_oi[pool_oi$tophit == dog_id,]
  if (nrow(mismatch_hit) > 0){
    print(mismatch_hit)
  }
}

#swaps in pid 10560
pool_oi <- mismatches[mismatches$pid == "PID_10560"
           ,]
for (dog_id in pool_oi$dog_id){
  mismatch_hit <- pool_oi[pool_oi$tophit == dog_id,]
  if (nrow(mismatch_hit) > 0){
    print(mismatch_hit)
  }
}

# metaData$lid_pid[metaData$lid_pid == "LID_1071346_PID_10560" & metaData$lid == "LID_1071346"] <- "_LID_1071349_PID_10560"

#swaps in pid 10562
pool_oi <- mismatches[mismatches$pid == "PID_10562"
           ,]
for (dog_id in pool_oi$dog_id){
  mismatch_hit <- pool_oi[pool_oi$tophit == dog_id,]
  if (nrow(mismatch_hit) > 0){
    print(mismatch_hit)
  }
}

write.csv(x = mismatches, file = "metadata_samples/mismatches.csv")

#remove all the mismatches that we were not able to identify
metaData <- metaData[!metaData$lid_pid %in% rownames(mismatches),]
rownames(metaData) <- metaData$lid_pid <- gsub("_LID" , "LID", metaData$lid_pid)

###########################################################################################################################
###########################################################################################################################
#####################    write this shit out    ###########################################################################
###########################################################################################################################
###########################################################################################################################

write.csv(file = "metadata_samples/metaData_swaps_lidpid_dog_id_tophit.csv", x = metaData_swaps)

metaData <- metaData[!metaData$lid_pid %in% metaData_swaps$lid_pid,]

write_rds("metadata_samples/ALL-DAP-metaData-240124.rds", x = metaData)

p1 <- metaData_[metaData_$Cohort == "precision_1",]
for (dog_id in p1$dog_id[duplicated(p1$dog_id)]){
  min_reads <- min(p1$reads[p1$dog_id==dog_id])
  p1 <- p1[!(p1$dog_id == dog_id & p1$reads == min_reads),]
}
dim(p1)
write_rds("metadata_samples/P1-DAP-metaData-240510.rds", x = p1)

p2 <- metaData_[metaData_$Cohort == "precision_2",]
for (dog_id in p2$dog_id[duplicated(p2$dog_id)]){
  min_reads <- min(p2$reads[p2$dog_id==dog_id])
  p2 <- p2[!(p2$dog_id == dog_id & p2$reads == min_reads),]
}
dim(p2)
write_rds("metadata_samples/P2-DAP-metaData-240510.rds", x = p2)

p3 <- metaData_[metaData_$Cohort == "precision_3",]
for (dog_id in p3$dog_id[duplicated(p3$dog_id)]){
  min_reads <- min(p3$reads[p3$dog_id==dog_id])
  p3 <- p3[!(p3$dog_id == dog_id & p3$reads == min_reads),]
}
dim(p3)
write_rds("metadata_samples/P3-DAP-metaData-240510.rds", x = p3)

triad <- metaData[!grepl("precision", metaData$Cohort),]
