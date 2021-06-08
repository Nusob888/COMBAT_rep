## Script to replicate generation of:
# contig_qc metrics
##Author: Bo Sun
##Lab: Bashford Rogers
##Institute: University of Oxford


setwd("/well/combat/users/vkh192/repertoire/")
.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/well/combat/users/vkh192/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos = 'http://cran.ma.imperial.ac.uk/')
library(Biostrings)
library(tidyverse)
library(parallel)
library(RSQLite)

#get cores
numCores = detectCores()

#Working directory library
wd <- list()

#Commonly used paths in my working directory
wd$output <- "/well/combat/users/vkh192/bcr_release_1/"
wd$R <- "/well/combat/users/vkh192/repertoire_code/R/Functions/"

#Get Santiagos sequencing metadata ----
repmeta <- readr::read_tsv("/well/combat/projects/repertoire/singlecell/samples.metadata")
repmeta$`Sample Name` <- gsub(".*_g","g",repmeta$`Sample Name`)
repmeta$`Sample Name` <- gsub("_.*","",repmeta$`Sample Name`)
repmeta <- repmeta %>% 
  select(`Sequencing ID`, `Sample Name`, `Flowcell ID`, `Lane Number`)
colnames(repmeta) <- c("seq_ID", "gplex", "flowcell", "lane_number")

#read in metadata from clustering and clinicalteams----
meta <- read.csv("/well/combat/projects/repertoire/shared/GEX_full_metadata.csv") %>% select(-X)

#Get meta data from repertoire SQLite ----
Bconn <- dbConnect(RSQLite::SQLite(), "/well/combat/projects/repertoire/imgtout/BCR/BCRdb/BCRIMGT.db") 
BCR10X <- RSQLite::dbGetQuery(Bconn, "SELECT `barcode`, `contig_id`, `is_cell`, `high_confidence`, `length`, `chain`, 
 `c_gene`, `full_length`, `reads`, `umis`, `seq_ID`, `gplex`, `group` FROM `10X_filtered_contig_annotations` WHERE `group` LIKE '%dehop'")
dbDisconnect(Bconn)

#Process BCRIMGT output with BcellR (unpublished package) ----
dir <- "/well/combat/projects/repertoire/imgtout/BCR/filtered_dehop/"
source("/well/combat/users/vkh192/repertoire_code/R/Functions/Parse_IMGT.R")
BCRIMGT <- GetIMGT(dir)
IMGT <- rbind(BCRIMGT$Heavy, BCRIMGT$Light) %>% 
  mutate(contig_id = paste(barcode, contig, sep="-"))  %>% 
  mutate(barcode = paste(barcode, "-1", sep = "")) %>% 
  select(-contig) %>%
  dplyr::mutate(seq_ID = sample)

#Merge meta ----
#merge repmeta
IMGT <- merge(IMGT, repmeta, by = c("seq_ID"), all.x=TRUE)

#merge 10X dataset
IMGT <- merge(IMGT, BCR10X, by=c("barcode", "contig_id", "gplex", "seq_ID"), all.x=TRUE)

#join metadata
IMGT$barcode_seq <- paste(IMGT$barcode, IMGT$gplex, sep = "-")
IMGT <- IMGT %>% 
  left_join(meta, by = c("barcode_seq"="barcode_seq", "barcode", "gplex"))

#filter out contigs not assigned to a COMBAT ID
IMGT <-IMGT %>% 
  filter(!is.na(COMBAT_ID))

#Assign incomplete VDJs to non-complete VDJs. One can opt to include these at the sacrifice of some inconsistency in how somatic hypermutations are calculated. However, for COMBAT repertoire analysis, we have opted for the most stringent threshold
IMGT <- IMGT %>% 
  mutate(vdj_complete= ifelse(v_partial_missing != 0, FALSE, ifelse(j_partial_missing != 0, FALSE, TRUE)))

#Parse chains. Please see script parsechains.R for explanation of procedure ----
source("/well/combat/users/vkh192/repertoire_code/R/Functions/parsechains.R")
chains <- Parsechains(IMGT, filter_productive=TRUE,format="IMGT", sampleid="gplex")
Parsed <- chains %>% 
  left_join(IMGT, by=c("barcode_seq", "contig_id", "locus", "umis"))
print(cutoff) #0.125
saveRDS(Parsed, "/well/combat/users/vkh192/repertoire/data/Parsed_metadata_v002.rds")
#Parsed <- readRDS("/well/combat/users/vkh192/repertoire/data/Parsed_metadata_v002.rds")

##Create table of contig QC metrics for annotation team----
#Select on heavy chain qc metrics as light chains are prone to ambient RNA contamination
contig_qc <- merge(IMGT, Parsed %>% select(barcode_seq, locus, contig_qc) %>% unique, by=c("barcode_seq", "locus"), all.x=TRUE) %>% 
  select(barcode_seq, locus, contig_id, contig_qc, vdj_complete)

write.csv(contig_qc, file=gzfile("/well/combat/users/vkh192/repertoire/data/bcr_release_v001_contig_qc.csv.gz"))
write.csv(contig_qc, file=gzfile("/well/combat/shared/B_chains_prelim/bcr_release_v001_contig_qc.csv.gz"))



#Call B cell clones ----
bcell_singlets <- Parsed %>% 
  filter(subset=="D") %>% 
  filter(vdj_complete=="TRUE") %>%
  filter(grepl("singleton|passed_qc",contig_qc))

#Seperate into Heavy and Light, cloneIDs will be called on Heavy chains
Heavy <- bcell_singlets %>% 
  filter(locus == "Heavy")
Light <- Parsed %>% 
  filter(subset=="D") %>% 
  filter(locus == "Light")

source("/well/combat/users/vkh192/repertoire_code/R/Functions/Get10Xclones_COMBAT.R")

clones <- Get10XclonesCOMBAT(Heavy)

saveRDS(clones, "/well/combat/users/vkh192/repertoire/data/new_clones_v2.rds")
clones <- readRDS("/well/combat/users/vkh192/repertoire/data/new_clones_v2.rds")
oldclones <- readRDS("/well/combat/users/vkh192/repertoire/data/new_clones.rds")

clones$Heavy %>% select(clone_per_replicate, clone_per_baseID)

#Merge Heavy and light
HC <- clones$Heavy[,c(5:9, 11:14, 16:58,189:194)]
colnames(HC) <- paste(colnames(HC), "HC", sep="_")
HC <- HC %>% rename(barcode_seq = barcode_seq_HC)
write.csv(HC, file=gzfile("/well/combat/users/vkh192/repertoire/data/bcr_release_v001_HC_table.csv.gz"))
write.csv(HC, file=gzfile("/well/combat/shared/B_chains_prelim/bcr_release_v001_HC_table.csv.gz"))

#48166 nrow
LC <- Light[,c(1:5, 7, 13:55, 186)]
colnames(LC) <- paste(colnames(LC), "LC", sep="_")
LC <- LC %>% rename(barcode_seq = barcode_seq_LC)
write.csv(LC, file=gzfile("/well/combat/users/vkh192/repertoire/data/bcr_release_v001_LC_table.csv.gz"))
write.csv(LC, file=gzfile("/well/combat/shared/B_chains_prelim/bcr_release_v001_LC_table.csv.gz"))


#check paired numbers, 48025 HCs paired with LC
HC$barcode_seq %in% (LC$barcode_seq %>% unique) %>% summary
LC$barcode_seq %>% unique %in% HC$barcode_seq %>% summary

#merge HC and LC
mergeanno <- merge(HC,LC, by=c("barcode_seq"), all.x=TRUE)

#sanity check numbers
mergeanno %>% filter(umi_ranks_LC ==2) %>% nrow #4517
mergeanno %>% filter(is.na(umi_ranks_LC)) %>% nrow #141
mergeanno %>% filter(umi_ranks_HC ==1) %>% filter(umi_ranks_LC ==1) %>% nrow #48025

#only keep top ranked LC umis for release to avoid duplication
mergeanno_rank1s<- mergeanno %>% filter(is.na(umi_ranks_LC) | umi_ranks_LC != 2) %>% as.tibble()
mergeanno_rank1s %>% nrow
write.csv(mergeanno_rank1s, file=gzfile("/well/combat/users/vkh192/repertoire/data/bcr_release_v001_chainmerge_table.csv.gz"))
write.csv(mergeanno_rank1s, file=gzfile("/well/combat/shared/B_chains_prelim/bcr_release_v001_chainmerge_table.csv.gz"))

#Additional file generated for annotation team to add umis to all associated barcodes
#need to add to readme
all_HC_cells <- Parsed %>% 
  filter(locus=="Heavy") %>%
  select(barcode_seq, contig_id, umis, pct_immunoglobin, umi_ranks) %>% 
  rename(umis_HC = umis)

all_LC_cells <- Parsed %>% 
  filter(locus=="Light") %>%
  select(barcode_seq, contig_id, umis, pct_immunoglobin, umi_ranks) %>% 
  rename(umis_LC = umis)

write.csv(all_HC_cells, file=gzfile("/well/combat/shared/B_chains_prelim/bcr_release_v001_umis_HC.csv.gz"))
write.csv(all_LC_cells, file=gzfile("/well/combat/shared/B_chains_prelim/bcr_release_v001_umis_LC.csv.gz"))







