##This script is to replicate the concatenation and conversion of "dehop" all_contig.fasta repertoire files from COMBAT study to IMGT compatible fasta files. 
#Author: Bo Sun
#Group: Bashford-Rogers
#Date: 2020-06-29

require(Biostrings)
require(tidyverse)
require(parallel)

#get cores
numCores = detectCores()

#Working directory library
wd <- list()

#Commonly used paths in my working directory
wd$data   <- "/well/combat/users/vkh192/repertoire/data/"
wd$output <- "/well/combat/users/vkh192/repertoire/out/"
wd$R <- "/well/combat/users/vkh192/repertoire/R/"

#Get meta
repmeta <- readr::read_tsv("/well/combat/projects/repertoire/singlecell/samples.metadata")
repmeta$gplex <- gsub(".*_g","g",repmeta$`Sample Name`)
repmeta$gplex<- gsub("_.*","",repmeta$gplex)

#Get directories
dir_list <- list.files(path = "/well/combat/projects/repertoire/singlecell/dehop/10X-vdj/", pattern = "\\.fasta$", full.names = TRUE, recursive = TRUE)
dir_list <- tibble(dir_list = dir_list)
all_contig_dir<- dir_list %>% filter(grepl("filtered_contig", dir_list))
head(all_contig_dir$dir_list)

#Create repdir file that will be used as a hub for all metadata
repdir <- list.files(path = "/well/combat/projects/repertoire/singlecell/dehop/10X-vdj/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
repdir <- tibble(annotation_dir = repdir)
repdir <- repdir %>% filter(grepl("filtered_contig_annotations", annotation_dir))
repdir <- repdir %>% mutate(seq_ID = gsub(".*//", "", annotation_dir)) %>% mutate(seq_ID = gsub("\\/.*", "", seq_ID))

#load all contig annotations to create rds file that we can query for isotype
repdir$all_contig_annotations <- lapply(repdir$annotation_dir, function(x){data.table::fread(x)})
repdir$numseqs<- sapply(repdir$all_contig_annotations, function(x){as.numeric(nrow(x))})
repdir$celltype <- sapply(repdir$all_contig_annotations, function(x){x %>% select(chain) %>% unique})
repdir <- repdir %>% left_join(repmeta, by = c("seq_ID" = "Sequencing ID"))
repdir<- repdir %>% rename(sample_name = `Sample Name`)
repdir<- repdir %>% mutate(celltype = ifelse(grepl("TCR", sample_name), "TCR", ifelse(grepl("BCR", sample_name),"BCR","unknown")))
repdir$fasta_dir <- all_contig_dir$dir_list

#Get DNAstringlist
fasta<- mclapply(repdir$fasta_dir, function(x){readDNAStringSet(x)}, mc.cores=numCores)

#Add suffix of seq ID
names(fasta) <- repdir$seq_ID

for (i in 1:length(fasta)) {
  names(fasta[[i]]) <- paste(names(fasta[[i]]), names(fasta[i]), sep = "-")
}

#sanity check
head(names(fasta[[1]]))

#separate by celltype
names(fasta) <- repdir$celltype
TCR <- fasta[names(fasta) == "TCR"]
BCR <- fasta[names(fasta) == "BCR"]

#remove names prior to concat
names(TCR) <- ""
names(BCR) <- ""

mergeTCR <- do.call(c, TCR)
mergeBCR <- do.call(c, BCR)

#Check lengths to decide need for chunking
length(mergeTCR)
length(mergeBCR)

#Summarise width metrics of seqs
summary(width(mergeTCR))
summary(width(mergeTCR))

#Function to split DNAstringsets into 1E06 chunks
chunkDNAstringset <- function(stringset){
  chunk <- 1000000
  n <- length(stringset)
  r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
  split <- split(stringset,r)
  return(split)
}

#Get chunks
TCRsplit <- chunkDNAstringset(mergeTCR)
BCRsplit <- chunkDNAstringset(mergeBCR)

#sanity check
TCRsplit[[1]]
BCRsplit[[5]]

#write fasta files

for (i in 1:length(TCRsplit)){
  filename <- paste(names(TCRsplit), "fasta", sep = "_COMBAT_TCR_all_contig_dehop.")
  writeXStringSet(TCRsplit[[i]], paste0(wd$output, filename[i]), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
}

for (i in 1:length(BCRsplit)){
  filename <- paste(names(BCRsplit), "fasta", sep = "_COMBAT_BCR_all_contig_dehop.")
  writeXStringSet(BCRsplit[[i]], paste0(wd$output, filename[i]), append=FALSE,
                  compress=FALSE, compression_level=NA, format="fasta")
}
