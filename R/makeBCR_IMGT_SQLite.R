##Script to replicate concatenation of BCR IMGT output for cellranger fasta alignments

#Author: Bo Sun
#Group: Bashford-Rogers
#Date: 2020-07-08

library(tidyverse)
library(parallel)
numCores <- detectCores()

##Prep directories----
BCRdir <- tibble::as.tibble(list.files(path = "/well/combat/projects/repertoire/imgtout/BCR/", pattern = "\\.txt$", full.names = TRUE, recursive = TRUE))
colnames(BCRdir) <- c("directory")
BCRdir <- BCRdir %>% 
  filter(!grepl("README.txt", directory)) %>% 
  mutate(names = gsub(".*/", "", directory)) %>% mutate(names = gsub(".txt", "", names)) %>% 
  mutate(group = gsub(".*filtered", "filtered", directory)) %>% mutate(group = gsub(".txt", "", group)) %>% 
  mutate(group = gsub("/.*", "", group)) %>%
  separate(group, c("cellranger_out", "index_corr"), "_") %>%
  mutate(index_corr = gsub("2","", index_corr))

#read in IMGT txt files as a column list of tables
BCRdir$tables <- mclapply(BCRdir$directory, function(x){readr::read_tsv(x)}, mc.cores = numCores)
BCRdir$tables <- mapply(cbind, BCRdir$tables, "names"=BCRdir$names, SIMPLIFY=F)
BCRdir$tables <- mapply(cbind, BCRdir$tables, "group"=BCRdir$index_corr, SIMPLIFY=F)

#Check names of correction assignments
BCRdir %>% select(index_corr) %>% unique

##Generate tibble of IMGTouts per index-correction method ----

#Parse default (non-index corrected IMGT outs)
defaultlist <- BCRdir %>% 
  filter(grepl("default", index_corr)) %>% 
  group_by(names) %>% 
  arrange(desc(names)) %>% 
  group_split() %>% 
  lapply(., function(x){data.table::rbindlist(x$tables)})

default<- list()

for (i in 1:length(defaultlist)){
  default[[paste0(defaultlist[[i]]$names %>% unique)]]<- defaultlist[[i]]
}

##Parse dehop (index corrected by lane IMGT outs)
dehoplist <- BCRdir %>% 
  filter(grepl("dehop$", index_corr)) %>% 
  group_by(names) %>% 
  arrange(desc(names)) %>% 
  group_split() %>% 
  lapply(., function(x){data.table::rbindlist(x$tables)})

dehop<- list()

for (i in 1:length(dehoplist)){
  dehop[[paste0(dehoplist[[i]]$names %>% unique)]]<- dehoplist[[i]]
}

#Remove BCRdir to save virtual memory

##Parse dehopsep (index corrected by lane IMGT outs)
dehopseplist <- BCRdir %>% 
  filter(grepl("dehopsep$", index_corr)) %>% 
  group_by(names) %>% 
  arrange(desc(names)) %>% 
  group_split() %>% 
  lapply(., function(x){data.table::rbindlist(x$tables)})

dehopsep<- list()

for (i in 1:length(dehopseplist)){
  dehopsep[[paste0(dehopseplist[[i]]$names %>% unique)]]<- dehopseplist[[i]]
}

##Merge to tibble
IMGTout<- tibble(IMGT_file = names(default),default=default, dehop=dehop, dehopsep=dehopsep)
IMGTout <- t(IMGTout)
colnames(IMGTout) <- IMGTout[1,]
IMGTout <- IMGTout[-1,]

##Create SQLitedb

library(RSQLite)
conn <- dbConnect(RSQLite::SQLite(), "/well/combat/projects/repertoire/imgtout/BCR/BCRdb/BCRIMGT.db")

lapply(IMGTout[,c("10_V-REGION-mutation-hotspots")], function(x){ dbWriteTable(conn,"filtered_V_REGION_mutation_hotspots", x, append = TRUE)})
lapply(IMGTout[,c("1_Summary")], function(x){ dbWriteTable(conn,"filtered_Summary", x, append = TRUE)})
lapply(IMGTout[,c("2_IMGT-gapped-nt-sequences")], function(x){ dbWriteTable(conn,"filtered_IMGT_gapped_nt_sequences", x, append = TRUE)})
lapply(IMGTout[,c("3_Nt-sequences")], function(x){ dbWriteTable(conn,"filtered_Nt_sequences", x, append = TRUE)})
lapply(IMGTout[,c("4_IMGT-gapped-AA-sequences")], function(x){ dbWriteTable(conn,"filtered_IMGT_gapped_AA_sequences", x, append = TRUE)})
lapply(IMGTout[,c("5_AA-sequences")], function(x){ dbWriteTable(conn,"filtered_AA_sequences", x, append = TRUE)})
lapply(IMGTout[,c("6_Junction")], function(x){ dbWriteTable(conn,"filtered_Junction", x, append = TRUE)})
lapply(IMGTout[,c("7_V-REGION-mutation-and-AA-change-table")], function(x){ dbWriteTable(conn,"filtered_V_REGION_mutation_AA_change", x, append = TRUE)})
lapply(IMGTout[,c("8_V-REGION-nt-mutation-statistics")], function(x){ dbWriteTable(conn,"filtered_V_REGION_nt_mutationstats", x, append = TRUE)})
lapply(IMGTout[,c("9_V-REGION-AA-change-statistics")], function(x){ dbWriteTable(conn,"filtered_V_REGION_AA_mutationstats", x, append = TRUE)})














