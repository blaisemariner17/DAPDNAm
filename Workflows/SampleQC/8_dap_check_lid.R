library(tidyverse)
library(parallel)
library(data.table)

dap_lid_pid_list <- list.files("/scratch/nsnyderm/dap_rrbs/check_out")

genotype_check <- function(x){
  a=read_delim(paste0("/scratch/nsnyderm/dap_rrbs/check_out/",
                      x),col_names = T,delim = " ")
  a= a%>%mutate(p_homhet=perc_het_consistent*perc_hom_consistent) %>%
    arrange(desc(p_homhet)) %>% slice(1:3) %>%
    summarize(hit1_score=p_homhet[1], hit2_score=p_homhet[2],
              hit3_score=p_homhet[3],
              tophit=SampleID[1], secondhit=SampleID[2],
              thirdhit=SampleID[3])
  a$lid=x
  return(a)}

# for (i in 1:length(dap_lid_pid_list)) {
#   genotype_check(dap_lid_pid_list[[i]])
# }

dap_rrbs_gtcheck <- rbindlist(
  mclapply(dap_lid_pid_list,
                   genotype_check,
                   mc.cores=4)
                   )

View(dap_rrbs_gtcheck)
colnames(dap_rrbs_gtcheck)[7]="lid_pid"
str_split_i(dap_rrbs_gtcheck$lid_pid,"_PID_",1)->dap_rrbs_gtcheck$lid

write_rds(dap_rrbs_gtcheck, file = "metadata_samples/dap_check_lid_res.rds")
