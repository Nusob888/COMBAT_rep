library(tidyverse)
dir <- list.files("/well/combat/publicDatasets/initial_deposition/CBD-KEY-CITESEQ-VDJ-B", full.names = TRUE)
dirnames <- list.files("/well/combat/publicDatasets/initial_deposition/CBD-KEY-CITESEQ-VDJ-B", full.names = FALSE)
dirnames <- gsub(".csv.gz", "",dirnames[-1])

df1 <- data.table::fread(dir[2])
df2 <- data.table::fread(dir[3])
df3 <- data.table::fread(dir[4])
df4 <- data.table::fread(dir[5])
df5 <- data.table::fread(dir[6])
df6 <- data.table::fread(dir[7])

comb <- list(df1, df2, df3, df4, df5, df6)
names(comb) <- dirnames

colnames(comb[[paste0(dirnames[2])]])

comb[[5]] %>% colnames

##define terms of colnames of sensitive data for restricted release
grepitems <- c("V1",
               "v_alignment", 
               "sequence_alignment_", 
               "d_j_region_",
               "v_region_", 
               "j_region_",
               "junction_nt_", 
               "sequence_", 
               "v_alignment_aa", 
               "sequence_alignment_aa",
               "vdj_complete"
               )
paste0(grepitems, collapse = "|")
test <- comb[[1]] %>% as_tibble()

comb2 <- lapply(comb, function(x){
  x <- as_tibble(x)
  x<- x[,!grepl(paste0(grepitems, collapse = "|"), names(x))]
  return(x)
})

names(comb2)
comb2[[6]]

lapply(names(comb2), function(x){
  write.csv(comb2[[x]], file=gzfile(paste0("/well/combat/publicDatasets/initial_deposition/CBD-KEY-CITESEQ-VDJ-B/",x,".csv.gz")))
})

#mv /well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-B/SEQUENCE_FILES /well/combat/publicDatasets/initial_deposition/CBD-RAW-REPERTOIRE-B-FASTA/
#mv /well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-T/SEQUENCE_FILES /well/combat/publicDatasets/initial_deposition/CBD-RAW-REPERTOIRE-T-FASTA/
  

