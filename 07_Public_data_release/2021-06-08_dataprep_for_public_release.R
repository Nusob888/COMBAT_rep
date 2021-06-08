##script to prepare public dataset release.
##Removal of clinical metadata
##Splitting of B and T cell data into respective directories
##Author: Bo Sun
##Lab: Bashford-Rogers
##Institute: University of Oxford

##Transfer of annotated files
dir <- list.files("/well/combat/datasets/CBD-BIR-00001/")

df1 <- data.table::fread(paste0("/well/combat/datasets/CBD-BIR-00001/", dir[2]))
df2 <- data.table::fread(paste0("/well/combat/datasets/CBD-BIR-00001/", dir[3]))
df3 <- data.table::fread(paste0("/well/combat/datasets/CBD-BIR-00001/", dir[4]))
df4 <- data.table::fread(paste0("/well/combat/datasets/CBD-BIR-00001/", dir[5]))

colnames(df1)
colnames(df2)
colnames(df3)
colnames(df4)
df1 <- df1[,-c(4:15)]
df2 <- df2[,-c(4:15)] 
df3 <- df3[,-c(4:15)] 
df4 <- df4[,-c(4:15)] 

colnames(df1)
colnames(df2)
colnames(df3)
colnames(df4)

listdf <- list(df1, df2, df3, df4)
names(listdf) <- gsub(".txt","", dir[2:5])

names(listdf[1:2])

lapply(names(listdf[1:2]), function(x){
  readr::write_csv(listdf[[x]], paste0("/well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-B/", x, ".csv"))
})

lapply(names(listdf[3:4]), function(x){
  readr::write_csv(listdf[[x]], paste0("/well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-T/", x, ".csv"))
})

##Transfer of IMGT files
dir2 <- list.files("/well/combat/datasets/CBD-BIR-00001/IMGT_SPLIT", full.names = TRUE)

BCRdir <- dir2[grepl("BCR",dir2)]
TCRdir <- dir2[grepl("TCR",dir2)]
BCRdir[1]
length(BCRdir)
length(TCRdir)

lapply(BCRdir, function(x){
  system(paste0("cp -p ",x,  " /well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-B/IMGT_SPLIT/"))
})
lapply(TCRdir, function(x){
  system(paste0("cp -p ",x,  " /well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-T/IMGT_SPLIT/"))
})

##Transfer of sequencing files
dir3 <- list.files("/well/combat/datasets/CBD-BIR-00001/SEQUENCE_FILES", full.names = TRUE)

BCRdir2 <- dir3[grepl("BCR",dir3)]
TCRdir2 <- dir3[grepl("TCR",dir3)]
BCRdir2[1]
length(BCRdir2)
length(TCRdir2)

lapply(BCRdir2, function(x){
  system(paste0("cp -p ",x,  " /well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-B/SEQUENCE_FILES/"))
})
lapply(TCRdir2, function(x){
  system(paste0("cp -p ",x,  " /well/combat/publicDatasets/initial_deposition/CBD-KEY-REPERTOIRE-T/SEQUENCE_FILES/"))
})
