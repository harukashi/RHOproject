setwd("./MB_SAMSI")
getwd()

# CRAN packages
packages <- c("tidyr","dplyr")
install.packages(packages, dep=TRUE, repos = "https://cloud.r-project.org/")

file.names=list.files(pattern = "*.rds")
list2env(lapply(setNames(file.names, make.names(gsub("*.rds$", "", gsub("^MB_SAMSI_", "", file.names)))), readRDS), envir = .GlobalEnv)
rm(file.names)

master=merge(subject_data,sample_data_stacked,"StudyID")
master=merge(master,rbind(dust_OTU,nasal_OTU),"sampleid")

saveRDS(master,"MB_SAMSI_master.rds")

master=readRDS("MB_SAMSI_master.rds")
