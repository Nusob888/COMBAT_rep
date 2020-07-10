##Script to replicate exploration of TCR BCR filtered outputs

#Author: Bo Sun
#Group: Bashford-Rogers
#Date: 2020-07-09

library(tidyverse)
library(RSQLite)
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

#load metadata
meta <- data.table::fread("/well/combat/projects/repertoire/shared/GEX_prelim_metadata.csv")
head(meta)

#Load IMGT ouput of choice. Summary usually has most of the interesting bits

##Prep TCR summary file ----
#make connection to database
conn <- dbConnect(RSQLite::SQLite(), "/well/combat/projects/repertoire/imgtout/TCR/TCRdb/TCRIMGT.db") 
#show list of available tables
dbListTables(conn)
#Extract 'dehop' IMGToutput from the Summary.txt output from IMGT
TCRsummary <- RSQLite::dbGetQuery(conn, "SELECT * FROM filtered_Summary WHERE `group` = 'dehopsep'")
head(TCRsummary)
#Separate Sequence ID into barcodes and seq_ID
TCRsummary <- data.table::as.data.table(TCRsummary %>% 
                                          separate(`Sequence ID`, c("barcode", "contig_id", "seq_ID"), "-") %>% 
                                          mutate(barcode = paste0(barcode, "-1", sep="")) %>% 
                                          mutate(contig_id = paste0(barcode,gsub("^1", "", contig_id), sep="")))
head(TCRsummary)
#First merge Santiagos metadata to IMGT by seq_ID
TCRsummary <- merge(TCRsummary, repmeta, by = c("seq_ID"), all.x=TRUE)
head(TCRsummary)
#sanity check the values against the original gplex numbers, here both sum to 70 unique combinations of seq_ID and gplex
TCRsummary %>% select(seq_ID, gplex) %>% unique %>% nrow
repmeta %>% select(seq_ID, sample_name, gplex) %>% filter(grepl("TCR", sample_name)) %>% unique %>% nrow

nrow(TCRsummary)

dbDisconnect(conn)

##Prep BCR summary file ----
conn <- dbConnect(RSQLite::SQLite(), "/well/combat/projects/repertoire/imgtout/BCR/BCRdb/BCRIMGT.db") 
#show list of available tables
dbListTables(conn)
#Extract 'dehop' IMGToutput from the Summary.txt output from IMGT
BCRsummary <- RSQLite::dbGetQuery(conn, "SELECT * FROM filtered_Summary WHERE `group` = 'dehopsep'")
head(BCRsummary)
#Separate Sequence ID into barcodes and seq_ID
BCRsummary <- data.table::as.data.table(BCRsummary %>% 
                      separate(`Sequence ID`, c("barcode", "contig_id", "seq_ID"), "-") %>% 
                      mutate(barcode = paste0(barcode, "-1", sep="")) %>% 
                      mutate(contig_id = paste0(barcode,gsub("^1", "", contig_id), sep="")))
head(BCRsummary)
#First merge Santiagos metadata to IMGT by seq_ID
BCRsummary <- merge(BCRsummary, repmeta, by = c("seq_ID"), all.x=TRUE)
head(BCRsummary)
#sanity check the values against the original gplex numbers, here both sum to 70 unique combinations of seq_ID and gplex
BCRsummary %>% select(seq_ID, gplex) %>% unique %>% nrow
repmeta %>% select(seq_ID, sample_name, gplex) %>% filter(grepl("BCR", sample_name)) %>% unique %>% nrow

dbDisconnect(conn)
colnames(BCRsummary)



TCRsummary <- data.table::as.data.table(TCRsummary)
BCRsummary <- data.table::as.data.table(BCRsummary)
##add 10X annotations to TCR and BCR summaries for comparison----
#Get directories
dir_list <- list.files(path = "/well/combat/projects/repertoire/singlecell/dehop.sep/10X-vdj/", pattern = "\\.fasta$", full.names = TRUE, recursive = TRUE)
dir_list <- tibble(dir_list = dir_list)
all_contig_dir<- dir_list %>% filter(grepl("filtered_contig", dir_list))
head(all_contig_dir$dir_list)
length(all_contig_dir$dir_list)

#Create repdir file that will be used as a hub for all metadata
repdir <- list.files(path = "/well/combat/projects/repertoire/singlecell/dehop.sep/10X-vdj/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
repdir <- tibble(annotation_dir = repdir)
repdir <- repdir %>% filter(grepl("filtered_contig_annotations", annotation_dir))
repdir <- repdir %>% mutate(seq_ID = gsub(".*//", "", annotation_dir)) %>% mutate(seq_ID = gsub("\\/.*", "", seq_ID))

#load all contig annotations to create rds file that we can query for isotype
repdir$all_contig_annotations <- mclapply(repdir$annotation_dir, function(x){data.table::fread(x)}, mc.cores=numCores)
repdir$numseqs<- sapply(repdir$all_contig_annotations, function(x){as.numeric(nrow(x %>% filter(grepl("IGH|TRB|TRD|Multi", chain))))})
repdir$celltype <- sapply(repdir$all_contig_annotations, function(x){x %>% select(chain) %>% unique})
repdir <- repdir %>% left_join(repmeta, by = c("seq_ID"))
repdir<- repdir %>% mutate(celltype = ifelse(grepl("TCR", sample_name), "TCR", ifelse(grepl("BCR", sample_name),"BCR","unknown")))
repdir %>% select(annotation_dir)
repdir <- repdir %>% mutate(post_process = gsub(".*singlecell/", "", annotation_dir)) %>% mutate(post_process = gsub("\\/10X.*", "", post_process))
repdir %>% select(post_process) %>% unique
repdir <- repdir %>% mutate(pool = gsub("[0-9]", "", gplex))
repdir$all_contig_annotations  <- mapply(cbind, repdir$all_contig_annotations, "seq_ID"=repdir$seq_ID, SIMPLIFY=F)
repdir$all_contig_annotations  <- mapply(cbind, repdir$all_contig_annotations, "gplex"=repdir$gplex, SIMPLIFY=F)


##Merge repdir annotations with summary files ----

#Get seq_IDs and filter repdir

#GetBCR annotations
BCRseq_ID <- BCRsummary$seq_ID %>% unique()
summary(repdir$seq_ID %in% BCRseq_ID)
BCRrepdir <- repdir[c(repdir$seq_ID %in% BCRseq_ID),]
BCRanno <- data.table::rbindlist(BCRrepdir$all_contig_annotations)
BCRanno <- BCRanno %>% select(barcode,contig_id, seq_ID, productive, is_cell, high_confidence, length, chain, v_gene, j_gene, c_gene, full_length, cdr3, reads, umis, gplex)
#add chain to BCRsummary
BCRsummary <- BCRsummary %>% mutate(chain = ifelse(grepl("IGH",`V-GENE and allele`), "IGH", ifelse(grepl("IGK", `V-GENE and allele`), "IGK", "IGL"))) 
BCRsummary <- data.table::as.data.table(BCRsummary)
#merge annotation and summary
BCRsummary <- merge(BCRsummary, BCRanno, by=c("gplex", "contig_id"), all.x=TRUE)
nrow(BCRsummary)
#Next merge the meta to by barcode and gPlex
BCRsummary <- merge(BCRsummary %>% rename(barcode = barcode.x), meta, by = c("barcode.x" = "barcode", "gplex"), all.x=TRUE)
nrow(BCRsummary)

#GetTCR annotations
TCRseq_ID <- TCRsummary$seq_ID %>% unique()
summary(repdir$seq_ID %in% TCRseq_ID)
TCRrepdir <- repdir[c(repdir$seq_ID %in% TCRseq_ID),]
TCRanno <- data.table::rbindlist(TCRrepdir$all_contig_annotations)
TCRanno <- TCRanno %>% select(barcode,contig_id, seq_ID, productive, is_cell, high_confidence, length, chain, v_gene, j_gene, c_gene, full_length, cdr3, reads, umis, gplex)
#add chain to TCR summary
TCRsummary <- TCRsummary %>% mutate(chain = ifelse(grepl("TRB",`V-GENE and allele`), "TRB", ifelse(grepl("TRA", `V-GENE and allele`), "TRA",ifelse(grepl("TRD", `V-GENE and allele`),"TRD" ,ifelse(grepl("TRG", `V-GENE and allele`), "TRG", "unknown")))))
TCRsummary <- data.table::as.data.table(TCRsummary)

#merge annotation and summary
TCRsummary <- merge(TCRsummary, TCRanno, by=c("gplex", "contig_id"), all.x=TRUE)
nrow(TCRsummary)
#Next merge the meta to by barcode and gPlex
TCRsummary <- merge(TCRsummary %>% rename(barcode = barcode.x), meta, by = c("barcode.x" = "barcode", "gplex"), all.x=TRUE)
nrow(TCRsummary)

#----

#remove duplicated barcodes that exist in the meta and those removed by QC
dedup <- meta %>% filter(!grepl("DOUBLET|UNASSIGNED", demuxletV2)) %>% filter(!is.na(cell_type))
dedup <- dedup[!c(dedup%>% select(barcode, gplex) %>% duplicated),] %>% data.table::as.data.table()
dedup %>% select(barcode, gplex) %>% duplicated %>% summary
dedup <- dedup %>% mutate(barcode_plex = paste0(barcode, gplex, sep=""))

#extract BCRs that were matched and unmatched to release 1 GEX barcodes
BCRsummary <- BCRsummary %>%mutate(barcode_plex = paste0(barcode, gplex, sep=""))
Bcellmatched <- BCRsummary[c(test1$barcode_plex %in% dedup$barcode_plex),] 

#extract TCRs that were matched and unmatched to release 1 GEX barcodes
TCRsummary <- TCRsummary %>%mutate(barcode_plex = paste0(barcode, gplex, sep=""))
Tcellmatched <- TCRsummary[c(test1$barcode_plex %in% dedup$barcode_plex),]

##Check barcode overlap per gplex ----
nrow(Bcellmatched)
nrow(Tcellmatched)

Bcellmatched$barcode_plex %in% Tcellmatched$barcode_plex %>% summary
Tcellmatched$barcode_plex %in% Bcellmatched$barcode_plex %>% summary


TCRdup <- Tcellmatched[c(Tcellmatched$barcode_plex %in% Bcellmatched$barcode_plex), ]
BCRdup <- Bcellmatched[c(Bcellmatched$barcode_plex %in% Tcellmatched$barcode_plex), ]
nrow(TCRdup)
nrow(BCRdup)

mergedup <- merge(TCRdup, BCRdup, by =c("barcode_plex"), all=TRUE)
nrow(mergedup)
head(mergedup)

mergedup %>% select(v_gene.x, v_gene.y, umis.x, umis.y, reads.x, reads.y, cell_type.x, cell_type.y) %>% head
#explore chain missassignments----


##All meta is now added. NAs are down to a combination of things: barcodes that do not exist in the GEX data either by change or inability to demultiplex and also from cellQC filtering 