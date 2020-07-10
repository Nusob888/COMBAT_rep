##Script to generate SQLite database for 10Xvdj annotations
#Author: Bo Sun
#Group: Bashford-Rogers
#Date: 2020-07-10

library(tidyverse)
library(parallel)
numCores <- detectCores()

#Get Santiago's metadata----
repmeta <- readr::read_tsv("/well/combat/projects/repertoire/singlecell/samples.metadata")
colnames(repmeta) 
repmeta$gplex <- gsub(".*_g","g",repmeta$`Sample Name`)
repmeta$gplex <- gsub("_.*","",repmeta$gplex)
repmeta <- repmeta %>% select(`Sequencing ID`, gplex, `Sample Name`, `Flowcell ID`, `Lane Number`)
colnames(repmeta) <- c("seq_ID", "gplex", "sample_name", "flowcell", "lane_number")
head(repmeta)


##add 10X annotations to TCR and BCR summaries for comparison----
#Get directories
dir_list <- list.files(path = "/well/combat/projects/repertoire/singlecell/", pattern = "\\.fasta$", full.names = TRUE, recursive = TRUE)
dir_list <- tibble(dir_list = dir_list)
all_contig_dir<- dir_list %>% filter(grepl("filtered_contig", dir_list))
all_contig_dir<- all_contig_dir %>% filter(!grepl("10X-vdj.dnv", dir_list))
all_contig_dir<- all_contig_dir %>% filter(grepl("dehop|dehopsep|default", dir_list))


head(all_contig_dir$dir_list)
length(all_contig_dir$dir_list)

#Create repdir file that will be used as a hub for all metadata
repdir <- list.files(path = "/well/combat/projects/repertoire/singlecell/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
repdir <- tibble(annotation_dir = repdir)
repdir <- repdir %>% filter(grepl("filtered_contig_annotations", annotation_dir))
repdir<- repdir %>% filter(!grepl("10X-vdj.dnv", annotation_dir))
repdir<- repdir %>% filter(grepl("dehop|dehopsep|default", annotation_dir))
repdir <- repdir %>% mutate(seq_ID = gsub(".*vdj/", "", annotation_dir)) %>% mutate(seq_ID = gsub("\\/.*", "", seq_ID))

#load all contig annotations to create rds file that we can query for isotype
repdir$all_contig_annotations <- mclapply(repdir$annotation_dir, function(x){data.table::fread(x)}, mc.cores=numCores)
repdir$numseqs<- sapply(repdir$all_contig_annotations, function(x){as.numeric(nrow(x %>% filter(grepl("IGH|TRB|TRD|Multi", chain))))})
repdir <- repdir %>% left_join(repmeta, by = c("seq_ID"))
repdir<- repdir %>% mutate(celltype = ifelse(grepl("TCR", sample_name), "TCR", ifelse(grepl("BCR", sample_name),"BCR","unknown")))
repdir %>% select(annotation_dir)
repdir <- repdir %>% mutate(post_process = gsub(".*singlecell//", "", annotation_dir)) %>% mutate(post_process = gsub("\\/10X.*", "", post_process)) %>%
  mutate(post_process = gsub("\\.", "", post_process))
repdir %>% select(post_process) %>% unique
repdir <- repdir %>% mutate(pool = gsub("[0-9]", "", gplex))

repdir$all_contig_annotations  <- mapply(cbind, repdir$all_contig_annotations, "seq_ID" = repdir$seq_ID, SIMPLIFY=F)
repdir$all_contig_annotations  <- mapply(cbind, repdir$all_contig_annotations, "gplex" = repdir$gplex, SIMPLIFY=F)
repdir$all_contig_annotations  <- mapply(cbind, repdir$all_contig_annotations, "group" = repdir$post_process, SIMPLIFY=F)

TCR <- repdir %>% filter(grepl("TCR", celltype))
BCR <- repdir %>% filter(grepl("BCR", celltype))

##Create SQLitedb

library(RSQLite)
conn <- dbConnect(RSQLite::SQLite(), "/well/combat/projects/repertoire/imgtout/TCR/TCRdb/TCRIMGT.db")
lapply(TCR$all_contig_annotations, function(x){ RSQLite::dbWriteTable(conn,"10X_filtered_contig_annotations", x, append = TRUE)})
dbListTables(conn)
dbDisconnect(conn)

conn <- dbConnect(RSQLite::SQLite(), "/well/combat/projects/repertoire/imgtout/BCR/BCRdb/BCRIMGT.db")
lapply(BCR$all_contig_annotations, function(x){ dbWriteTable(conn,"10X_filtered_contig_annotations", x, append = TRUE)})
dbListTables(conn)
dbDisconnect(conn)
