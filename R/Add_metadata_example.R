##Script to replicate merging of TCR IMGT output with metadata

#Author: Bo Sun
#Group: Bashford-Rogers
#Date: 2020-07-08

library(tidyverse)
library(RSQLite)

#Get Santiago's metadata
repmeta <- readr::read_tsv("/well/combat/projects/repertoire/singlecell/samples.metadata")
colnames(repmeta) 
repmeta$gplex <- gsub(".*_g","g",repmeta$`Sample Name`)
repmeta$gplex <- gsub("_.*","",repmeta$gplex)
repmeta <- repmeta %>% select(`Sequencing ID`, gplex, `Sample Name`, `Flowcell ID`, `Lane Number`)
colnames(repmeta) <- c("seq_ID", "gplex", "sample_name", "flowcell", "lane_number")
head(repmeta)

#load metadata
meta <- data.table::fread("/well/combat/projects/repertoire/shared/GEX_prelim_metadata.csv")
head(meta)

#Load IMGT ouput of choice. Summary usually has most of the interesting bits

#make connection to database
conn <- dbConnect(RSQLite::SQLite(), "/well/combat/projects/repertoire/imgtout/TCR/TCRdb/TCRIMGT.db") 

#show list of available tables
dbListTables(conn)

#Extract 'dehop' IMGToutput from the Summary.txt output from IMGT
Summary <- RSQLite::dbGetQuery(conn, "SELECT * FROM filtered_Summary WHERE `group` = 'dehop'")
head(Summary)

#Separate Sequence ID into barcodes and seq_ID
Summary <- data.table::as.data.table(Summary %>% separate(`Sequence ID`, c("barcode", "contig", "seq_ID"), "-") %>% mutate(barcode = paste0(barcode, "-1", sep="")))
head(Summary)

#First merge Santiagos metadata to IMGT by seq_ID
Summary <- merge(Summary, repmeta, by = c("seq_ID"), all.x=TRUE)
head(Summary)

#sanity check the values against the original gplex numbers, here both sum to 70 unique combinations of seq_ID and gplex
Summary %>% select(seq_ID, gplex) %>% unique %>% nrow
repmeta %>% select(seq_ID, sample_name, gplex) %>% filter(grepl("TCR", sample_name)) %>% unique %>% nrow

#Next merge the meta to by barcode and gPlex
Summary <- merge(Summary, meta, by = c("barcode", "gplex"), all.x=TRUE)
head(Summary)

##All meta is now added. NAs are down to a combination of things: barcodes that do not exist in the GEX data either by change or inability to demultiplex and also from cellQC filtering 