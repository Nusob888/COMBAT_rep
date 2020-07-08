##Script to replicate 10Xvdj output comparison

#Author: Bo Sun
#Group: Bashford-Rogers
#Date: 2020-06-29

setwd("/well/combat/users/vkh192/repertoire/")
.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/well/combat/users/vkh192/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos = 'http://cran.ma.imperial.ac.uk/')
library(Biostrings)
library(tidyverse)
library(parallel)

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
colnames(repmeta) 

repmeta$gplex <- gsub(".*_g","g",repmeta$`Sample Name`)
repmeta$gplex<- gsub("_.*","",repmeta$gplex)

head(repmeta$gplex)

#Get directories
dir_list <- list.files(path = "/well/combat/projects/repertoire/singlecell/", pattern = "\\.fasta$", full.names = TRUE, recursive = TRUE)
dir_list <- tibble(dir_list = dir_list)
all_contig_dir<- dir_list %>% filter(!grepl("dnv", dir_list))
all_contig_dir<- dir_list %>% filter(grepl("filtered_contig", dir_list))
head(all_contig_dir$dir_list)
length(all_contig_dir$dir_list)

#Create repdir file that will be used as a hub for all metadata
repdir <- list.files(path = "/well/combat/projects/repertoire/singlecell/", pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
repdir <- tibble(annotation_dir = repdir)
repdir<- repdir %>% filter(!grepl("dnv", annotation_dir))
repdir <- repdir %>% filter(grepl("filtered_contig_annotations", annotation_dir))
repdir <- repdir %>% mutate(seq_ID = gsub(".*vdj/", "", annotation_dir)) %>% mutate(seq_ID = gsub("\\/.*", "", seq_ID))

#load all contig annotations to create rds file that we can query for isotype
repdir$all_contig_annotations <- mclapply(repdir$annotation_dir, function(x){data.table::fread(x)}, mc.cores=numCores)
repdir$numseqs<- sapply(repdir$all_contig_annotations, function(x){as.numeric(nrow(x %>% filter(grepl("IGH|TRB|TRD|Multi", chain))))})
repdir$celltype <- sapply(repdir$all_contig_annotations, function(x){x %>% select(chain) %>% unique})
repdir <- repdir %>% left_join(repmeta, by = c("seq_ID" = "Sequencing ID"))
repdir<- repdir %>% rename(sample_name = `Sample Name`)
repdir<- repdir %>% mutate(celltype = ifelse(grepl("TCR", sample_name), "TCR", ifelse(grepl("BCR", sample_name),"BCR","unknown")))

repdir %>% select(annotation_dir)
repdir <- repdir %>% mutate(post_process = gsub(".*//", "", annotation_dir)) %>% mutate(post_process = gsub("\\/10X.*", "", post_process))

repdir %>% select(post_process) %>% unique
repdir <- repdir %>% mutate(gpool = gsub("[0-9]", "", gplex))

pdf("/well/combat/users/vkh192/repertoire/out/dehopvsdefault.pdf", height = 10, width = 10)
repdir %>% select(`Lane Number`, `seq_ID`, `gplex`, `post_process`, `numseqs`, celltype) %>%
  ggplot(aes(gplex, numseqs, fill=post_process)) +
  geom_col(position="dodge")+
  facet_wrap(~celltype)+
  coord_flip()
dev.off()

pdf("/well/combat/users/vkh192/repertoire/out/dehopvsdefault_gplexpool.pdf", height = 10, width = 10)
repdir %>% select(`Lane Number`, `seq_ID`, `gplex`, `post_process`, `numseqs`, celltype, gpool) %>%
  ggplot(aes(gpool, numseqs, fill=post_process)) +
  geom_col(position="dodge")+
  facet_wrap(~celltype)+
  coord_flip()
dev.off()

saveRDS(repdir, paste0(wd$data, "repdir_allhopconds.rds"))

repdir <- readRDS(paste0(wd$data, "repdir_allhopconds.rds"))

#get meta -----

library(RSQLite)
conn <- RSQLite::dbConnect(SQLite(), "/well/combat/projects/preprocess/citeseq_initial/pipeline_celldb/170620/csvdb")
#to see what tables there are dbListTables(conn)
#to read use t <- dbReadTable(conn, "final")
meta <- dbReadTable(conn, "final")
head(meta)

#Get preprocessing team's preliminary clustering
cluster_barcodes<- readr::read_tsv("/well/combat/projects/citeseq/data/gex/lowdepth/1.0/baseline_analysis/harmony.seurat.dir/50_1_1_wilcox/cluster.dir/cluster_assignments.txt.gz")
colnames(cluster_barcodes) <- c("barcode_id", "cluster_id")
head(cluster_barcodes)

cluster_annotation <- readr::read_tsv("/well/combat/projects/preprocess/citeseq_initial/pipeline_celldb/300620/cluster_info.txt")
head(cluster_annotation)

#join dataframes 
cluster <- merge(cluster_barcodes, cluster_annotation, by = c("cluster_id"), all.x=TRUE)
head(cluster)

clusters <- merge(cluster, meta, by = c("barcode_id"), all.x=TRUE)
head(clusters)

clusters <- clusters %>% group_by(cell_type, pool) %>% add_tally(name="cell_nums")
clusters <- clusters %>% select(cell_type, pool, cell_nums) %>% mutate(celltype= ifelse(grepl("T-cell", cell_type), "TCR", ifelse(grepl("B-cell", cell_type), "BCR", "Other")))
clusters <- clusters %>% filter(grepl("TCR|BCR", celltype))

colnames(clusters) <- c("cell_anno", "gpool", "numseqs", "celltype")
clusters$post_process <- c("prelim_clustering")


clusters <- rbind((clusters %>% ungroup %>% select(gpool, post_process, numseqs, celltype)), (repdir %>% select(gpool, post_process, numseqs, celltype)))


pdf("/well/combat/users/vkh192/repertoire/out/prelim_cluster_num.pdf", height = 6, width = 10)
clusters %>%
  ungroup() %>% 
  group_by(gpool, post_process, celltype) %>%
  unique() %>%
  ggplot(aes(gpool, numseqs, fill= post_process))+
  geom_col(position="dodge")+
  facet_wrap(~celltype)+
  coord_flip()
dev.off()



#per COMBAT sample distribution 
cluster_sample <- merge(cluster_barcodes, cluster_annotation, by = c("cluster_id"), all.x=TRUE)
head(cluster_sample)

cluster_sample <- merge(cluster_sample, meta, by = c("barcode_id"), all.x=TRUE)
head(cluster_sample)
nrow(cluster_sample)

cluster_sample <- cluster_sample %>% group_by(cell_type, pool, COMBATID) %>% add_tally(name="cell_nums")
cluster_sample <- cluster_sample %>% select(cell_type, pool, cell_nums, COMBATID, source) %>% mutate(celltype= ifelse(grepl("T-cell", cell_type), "TCR", ifelse(grepl("B-cell", cell_type), "BCR", "Other")))
cluster_sample <- cluster_sample %>% filter(grepl("TCR|BCR", celltype))

pdf("/well/combat/users/vkh192/repertoire/out/prelim_cluster_num_id.pdf", height = 8, width = 8)
cluster_sample %>%
  ungroup() %>% 
  group_by(COMBATID) %>%
  unique %>%
  ggplot(aes(COMBATID, cell_nums, fill=COMBATID))+
  geom_col(position= "dodge")+
  facet_grid(vars(source), vars(cell_type), scales="free")+
  coord_flip()+
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust = 1), axis.text.y = element_blank())
dev.off()


##Analysis per individual ----

repdir %>% select(post_process) %>% unique

TCR <- repdir %>% filter(grepl("dehop.sep", post_process)) %>% filter(grepl("TCR", celltype))
BCR <- repdir %>% filter(grepl("dehop.sep", post_process)) %>% filter(grepl("BCR", celltype))

#Add per seq_ID annotations
BCR$all_contig_annotations  <- mapply(cbind, BCR$all_contig_annotations, "seq_ID"=BCR$seq_ID, SIMPLIFY=F)
BCR$all_contig_annotations  <- mapply(cbind, BCR$all_contig_annotations, "gplex"=BCR$gplex, SIMPLIFY=F)

BCRdf <- data.table::rbindlist(BCR$all_contig_annotations)

##NOTE clonotypes are inaccurate due to the way cellranger was run. need to recalculate
saveRDS(BCRdf, paste0(wd$data, "BCRdf.rds"))

#create list of VJ pairs (only as IMGT hasnt finished)
clus_meta <- merge(BCRdf, meta, by = c("barcode", "gplex"), all.x=TRUE)
head(clus_meta)
head(cluster)
BCRdf <- merge(clus_meta, cluster, by = c("barcode_id"), all.x=TRUE)

#finer celltypes
BCRdf %>% ungroup() %>% group_by(name) %>% add_tally() %>% select(cell_type, n) %>% unique %>% as.data.frame()

#add sequence number and cleanup
BCRdf$sequence_number <- c(1:nrow(BCRdf))
BCRdf <- BCRdf %>% filter(productive == "TRUE")
Heavy <- BCRdf %>% filter(chain == "IGH") %>% filter(!is.na(COMBATID))
Light <- BCRdf %>% filter(grepl("IGL|IGK", chain)) %>% filter(!is.na(COMBATID))

#Check non assigned BCRs

pdf("/well/combat/users/vkh192/repertoire/out/umis_baseID.pdf", height = 4, width = 5)
BCRdf %>% filter(chain == "IGH") %>% mutate(baseID= ifelse(grepl("COMBAT", baseID), "singlet", baseID)) %>% 
  ggplot(aes(umis, fill=baseID))+
  geom_density()+
  scale_x_log10()+
  facet_wrap(~baseID, nrow = 1)+
  coord_flip()
dev.off()

pdf("/well/combat/users/vkh192/repertoire/out/reads_baseID.pdf", height = 4, width = 5)
BCRdf %>% filter(chain == "IGH") %>% mutate(baseID= ifelse(grepl("COMBAT", baseID), "singlet", baseID)) %>% 
  ggplot(aes(reads, fill=baseID))+
  geom_density()+
  scale_x_log10()+
  facet_wrap(~baseID, nrow = 1)+
  coord_flip()
dev.off()

pdf("/well/combat/users/vkh192/repertoire/out/contig_length.pdf", height = 4, width = 5)
BCRdf %>% filter(chain == "IGH") %>% mutate(baseID= ifelse(grepl("COMBAT", baseID), "singlet", baseID)) %>% 
  ggplot(aes(length, fill=baseID))+
  geom_density()+
  scale_x_log10()+
  facet_wrap(~baseID, nrow = 1)+
  coord_flip()
dev.off()


#reassign clones
clusterlist<- Heavy %>% ungroup() %>% group_by(v_gene, j_gene) %>% group_split()

#function to calculate cdr3dist ----
stringdistclusters <- function(clusterlist){
  numCores <- detectCores()
  clusterlist <- lapply(clusterlist, function(x){
    x$length <- nchar(x$cdr3_nt)
    dist <- stringdist::stringdistmatrix(x$cdr3_nt, x$cdr3_nt, method="lv")
    rownames(dist) <- x$sequence_number
    colnames(dist) <- x$sequence_number
    dist[upper.tri(dist)] <- NA
    dist <- reshape2::melt(dist)
    dist <- filter(dist, !is.na(value))
    colnames(dist) <- c("sequence_number", "seq_to", "lv_dist")
    dist$sequence_number <- as.numeric(dist$sequence_number)
    dist$seq_to <- as.numeric(dist$seq_to)
    dist <- dist %>% left_join(x[,c("sequence_number", "length")], by = "sequence_number") %>% dplyr::rename("s1_length" = "length")
    dist <- dist %>% left_join(x[,c("sequence_number", "length")], by = c("seq_to"="sequence_number")) %>% dplyr::rename("s2_length" = "length")
    dist <- dist %>% mutate(norm_lvdist = (2*lv_dist)/(s1_length+s2_length+lv_dist))
    x <- dist %>% select(-s1_length, -s2_length, cdr3_nt)
  })
  edgelist <- data.table::rbindlist(clusterlist) 
  return(edgelist)
}


#Function to assign clusters ----
library(igraph)
#function to calculate levenshtein distances of cdr3s after VJ clustering
igraphclusters <- function(cdr3_dist){
  numCores <- detectCores()
  #create list of samples
  cdr3list <- split(cdr3_dist, cdr3_dist$COMBATID)
  igraphlist <- list()
  
  cdr3list <- mclapply(cdr3list, function(x){
    distgraph <- igraph::graph.data.frame(x[,1:2])
    links <- data.table::data.table(sequence_number=unique(unlist(x[,1:2])), clone= paste(igraph::clusters(distgraph)$membership, unique(unlist(x[,c("COMBATID")])), sep = "-"))
    x<- x%>% left_join(links, by= c("sequence_number"))
  }, mc.cores = numCores) 
  
  clones <- data.table::rbindlist(cdr3list) 
  return(clones)
}


#calculate cdr3dist ----
cdr3dist<- stringdistclusters(clusterlist)

#Calculate the density cutoff
d <- density(cdr3dist$norm_lvdist)
cutoff<- optimize(approxfun(d$x,d$y),interval=c(0,0.4))$minimum

#add 1 to lv_dist for network plotting; distances of 0, where sequences are identical will not form edges
cdr3dist <- cdr3dist %>% mutate(norm_lvdist = norm_lvdist+1)

#threshold clone calls
cdr3dist <- cdr3dist %>% 
  dplyr::mutate(norm_lvdist = ifelse(norm_lvdist <cutoff+1, norm_lvdist, 0)) %>% 
  dplyr::filter(!norm_lvdist == 0)

#merge back to input df
cdr3dist <- cdr3dist %>% left_join(Heavy, by = "sequence_number") 
cdr3dist$sequence_number <- as.character(cdr3dist$sequence_number)
cdr3dist$seq_to <- as.character(cdr3dist$seq_to)

#add clones
cdr3dist<- igraphclusters(cdr3dist)

edgelist <- cdr3dist[,c("sequence_number", "seq_to", "cdr3_nt")]
nam <- cdr3dist %>% ungroup() %>% group_by(sequence_number) %>% dplyr::slice(1) %>% ungroup()

#Get clone abunances
nam <- nam %>% 
  dplyr::group_by(COMBATID) %>% 
  dplyr::add_tally(name= "tot_samp") %>% 
  dplyr::ungroup() %>% 
  dplyr::group_by(clone, COMBATID) %>% 
  dplyr::add_tally(name = "clonal_abundance") %>% 
  dplyr::mutate(seq_freq = clonal_abundance /tot_samp ) %>% 
  dplyr::ungroup() %>%
  dplyr::select(-tot_samp)

out<- list(nam, edgelist)
names(out) <- c("Heavy", "edgelist")

out
#----

out$Heavy %>% select(cell_type)

#check how many BCRs were lost/assigned to other celltypes
library(RColorBrewer)
# Define colours
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

pdf("/well/combat/users/vkh192/repertoire/out/BCR_cellassignment.pdf", height = 5, width = 5)
out$Heavy %>% filter(!is.na(cell_type)) %>% group_by(cell_type, source) %>% add_tally() %>% select(cell_type, source, n) %>% unique %>%
  ggplot(aes(source, n, fill=cell_type))+
  geom_col(position="fill")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_fill_manual(values= mycolors)
dev.off()

library(gghighlight)
pdf("/well/combat/users/vkh192/repertoire/out/clonal_abundance.pdf", height = 5, width = 5)
out$Heavy %>% filter(!is.na(cell_type)) %>% filter(grepl("B-cell|plasma", cell_type)) %>% 
  ungroup() %>% group_by(source) %>% add_tally(name="total_seqnum") %>% ungroup() %>%
  group_by(clone, source) %>% add_tally() %>% select(clone, source, n, total_seqnum) %>% unique() %>%
  arrange(desc(n)) %>%
  mutate(percent = n/total_seqnum) %>%
  ggplot(aes(source, percent, fill=clone))+
  geom_col(position="stack")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  gghighlight(n> 40) +
  scale_fill_simpsons()
dev.off()


#shared clones
shared_clones <- out$Heavy %>% group_by(cdr3_nt) %>% add_tally(name="seqfreq") %>% 
  ungroup() %>% group_by(COMBATID, cdr3_nt) %>% add_tally(name="sample_seqfreq") %>% 
  ungroup() %>% filter(!(sample_seqfreq %in% seqfreq))

knitr::kable(shared_clones %>%
               select(cell_type, COMBATID, cdr3, umis, reads) %>% 
               group_by(cell_type, COMBATID, cdr3) %>% 
               add_tally() %>% 
               mutate(mean_umis = median(umis)) %>%
               mutate(mean_reads = median(reads)) %>% 
               select(cell_type, COMBATID, cdr3,n, mean_umis, mean_reads) %>%
               unique %>% 
               arrange(desc(cdr3)))

#Identify trends of doublets
pdf("/well/combat/users/vkh192/repertoire/out/umis_celltypehist.pdf", height = 4, width = 11)
out$Heavy %>% filter(!is.na(cell_type)) %>%
  ggplot(aes(umis, fill=cell_type))+
  geom_density()+
  scale_x_log10()+
  facet_wrap(~cell_type, nrow = 1)+
  coord_flip()
dev.off()

pdf("/well/combat/users/vkh192/repertoire/out/reads_celltypehist.pdf", height = 4, width = 11)
out$Heavy %>% filter(!is.na(cell_type)) %>%
  ggplot(aes(reads, fill=cell_type))+
  geom_density()+
  scale_x_log10()+
  facet_wrap(~cell_type, nrow = 1)+
  coord_flip()
dev.off()
  
#check alleles
out$Heavy %>% separate(v_gene, c("v_fam", "allele"),"\\*") %>% filter(grepl("IGHM", c_gene)) %>% ungroup() %>%
  group_by(COMBATID, v_fam, allele) %>% slice(1) %>% ungroup() %>% group_by(COMBATID, v_fam) %>% add_tally(name="num_alleles") %>%
  select(COMBATID, v_fam, allele, num_alleles) %>% arrange(desc(num_alleles))

pdf("/well/combat/users/vkh192/repertoire/out/isotype.pdf", height = 5, width = 5)
out$Heavy %>% filter(!is.na(cell_type)) %>% filter(grepl("B-cell|plasma", cell_type)) %>% 
  ungroup() %>% group_by(source) %>% add_tally(name="total_seqnum") %>% ungroup() %>% mutate(c_gene = gsub("\\*.*", "", c_gene)) %>%
  group_by(c_gene, source) %>% add_tally() %>% select(c_gene, source, n, total_seqnum) %>% unique() %>%
  arrange(desc(n)) %>%
  mutate(proportion = n/total_seqnum) %>%
  ggplot(aes(source, proportion, fill=c_gene))+
  geom_col(position="stack")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_fill_d3()
dev.off()


pdf("/well/combat/users/vkh192/repertoire/out/v_gene.pdf", height = 13, width = 13)
out$Heavy %>% filter(!is.na(cell_type)) %>% filter(grepl("B-cell|plasma", cell_type)) %>% 
  ungroup() %>% group_by(source) %>% add_tally(name="total_seqnum") %>% ungroup() %>% mutate(c_gene = gsub("\\*.*", "", c_gene)) %>%
  separate(v_gene, c("v_gene", "allele"),"\\*") %>%
  group_by(v_gene, source) %>% add_tally() %>% select(v_gene, source, n, total_seqnum) %>% unique() %>%
  arrange(desc(n)) %>%
  mutate(proportion = n/total_seqnum) %>% 
  ggplot(aes(v_gene, proportion, fill=source))+
  geom_col(position="dodge")+
  scale_fill_d3()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  coord_flip()+
  facet_wrap(~source, nrow=1)
dev.off()

