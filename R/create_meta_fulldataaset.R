## This script is to replicate the generation of the metadata csv for the COMBAT repertoire team
## Author: Bo Sun
## Date: 2020-08-11
## Group: Bashford-Rogers

setwd("/well/combat/users/vkh192/repertoire/")
.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/well/combat/users/vkh192/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos = 'http://cran.ma.imperial.ac.uk/')
library(tidyverse)

#Get metadata from preprocessing SQLite db

repmeta <- readr::read_tsv("/well/combat/projects/repertoire/singlecell/samples.metadata")
colnames(repmeta) 
repmeta$`Sample Name` <- gsub(".*_g","g",repmeta$`Sample Name`)
repmeta$`Sample Name` <- gsub("_.*","",repmeta$`Sample Name`)
repmeta <- repmeta %>% select(`Sequencing ID`, `Sample Name`, `Flowcell ID`, `Lane Number`)

head(repmeta)
colnames(repmeta) <- c("seq_ID", "gplex", "flowcell", "lane_number")

#check new directory sampleqc
celltable <- readr::read_tsv("/well/combat/datasets/CBD-10X-00004/full_dataset/cell.table.tsv.gz")
celltable %>% select(scRNASeq_sample_ID)
celltable <- celltable %>% mutate(barcode_seq = paste(barcode, sequencing_id, sep="-"))

#Read in cell annotations
cluster <- readr::read_tsv("/well/combat/datasets/CBD-10X-00004/full_dataset/cluster.assignments.tsv.gz")
annotations <- readr::read_tsv("/well/combat/datasets/CBD-10X-00004/full_dataset/cluster.annotations.tsv")

#merge tables
metatable <- celltable %>% left_join(cluster, by=c("barcode_seq"="barcode"))

metatable <- metatable %>% left_join(annotations, by=c("cluster_id"))

#read clinical data
clin <- readr::read_tsv("/well/combat/datasets/clinicalData/basicData/CBD-CLIN-00004/Final_COMBAT_basic_clinical_data_freeze_210720.txt")
clin %>% as.data.frame() %>% head

#merge clinical data 
metatable <- metatable %>% left_join(clin, by=c("scRNASeq_sample_ID"))

write.csv(metatable, "/well/combat/projects/repertoire/shared/GEX_full_metadata.csv")



