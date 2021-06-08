##Script to process data for handover to annotation team. 
##Author: Bo Sun
##Lab: Bashford Rogers
##Institute: University of Oxford

setwd("/well/combat/users/vkh192/repertoire/")
.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/well/combat/users/vkh192/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos = 'http://cran.ma.imperial.ac.uk/')

library(tidyverse)
library(parallel)
library(Hmisc)

numCores <- detectCores()
##Data loading and cleaning ----
clones <- readRDS("/well/combat/users/vkh192/repertoire/data/new_clones.rds")

Heavy <- clones$Heavy
Light <- readRDS("/well/combat/users/vkh192/repertoire/data/Light_chain.rds")
Heavy$Source <- factor(Heavy$Source, levels = c("HV", "Flu","Sepsis", "COVID_LDN","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

umap <- data.table::fread("/well/combat/datasets/CBD-10X-00004/cell_subsets/d_bcell_plasma/umap.tsv.gz")
head(umap)
Heavy <- Heavy %>% left_join(umap,by=c("barcode_seq"="barcode"))

#load new annotations
newanno <- readr::read_csv("/well/combat/users/vkh192/repertoire/data/newanno.csv") %>% select(-name)
Heavy <- Heavy %>% left_join(newanno, by=c("barcode_seq"="X1"))
Light <- Light %>% left_join(newanno, by=c("barcode_seq"="X1"))

head(Heavy %>% select(new_clusters))
##Get clonality measures----

#Clonal expansion index is the simpsons index for each unique VDJ identity with uniform subsampling to the same number of unique BCRs

IDlist<- Heavy %>% 
  group_by(COMBAT_ID_Time, new_clusters) %>%
  group_split()

length(IDlist)
IDlist[[1]] %>% select(COMBAT_ID_Time, Source,new_clusters)
library(vegan)
library(mosaic)
library(DescTools)

#For CEI, we must first group by barcodes and VDJ sequence per sample to collapse down to unique. Resampling will be done to the lowest sample size
Heavy %>% 
  select(COMBAT_ID_Time, new_clusters) %>%
  group_by(COMBAT_ID_Time, new_clusters) %>% 
  add_tally() %>%
  ungroup%>%
  select(COMBAT_ID_Time, new_clusters, n)%>% 
  unique %>%
  as.data.frame() %>%
  arrange(n) %>% head
#Here we see the smallest number of individual BCRs is 25, therefore subsampling must capture at least this many. 
#If we replace during subsampling, we can still estimate diversity indirectly vs. sampling without replacement which only captures the same diversity

getCEI <- function(input, B = 10000 , conf.level = 0.95) {
  data <- input[,c("sequence")]
  counts <- table(data)
  ID <- input %>% select(COMBAT_ID_Time, Age, Sex ,Source, new_clusters) %>% unique
  options(`mosaic:parallelMessage` = FALSE)
  gini.boot <- mosaic::do(B) * Gini(resample(counts, 20, replace=TRUE))
  names(gini.boot) <- "result"
  boot.mean <- mean(~result, data = gini.boot)
  boot.se <- sd(~result, data = gini.boot)
  boot.median <- median(~result, data = gini.boot)
  g <- as.numeric(1 - conf.level)
  boot.lci <- mosaic::qdata(~result, g/2,data = gini.boot)
  boot.uci <- mosaic::qdata(~result,1-g, data = gini.boot)
  boot.stats <- data.table::data.table(mean=boot.mean, se=boot.se, media=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
  boot.stats <- cbind(ID, boot.stats)
  bootdata <- data.table::data.table(ID, gini_boot=gini.boot)
  return(bootdata)
}

CEI <- mclapply(IDlist, function(x){getCEI(input=x)}, mc.cores=numCores)
CEI <- data.table::rbindlist(CEI)
CEImean <- aggregate(gini_boot.result ~ COMBAT_ID_Time+Source+new_clusters, data=CEI, mean) %>%
  rename(mean_cei=gini_boot.result)
CEImean %>%head

CEImean$Source <- factor(CEImean$Source, levels = c("HV", "Flu","Sepsis", "COVID_LDN","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

png("/well/combat/users/vkh192/repertoire/out/subset_cei.png",type=c("cairo"), height = 5, width = 8, units = 'in', res = 400)
CEImean %>%
  ggplot(aes(new_clusters, mean_cei, fill=Source, color=Source))+
  geom_boxplot(alpha=0.4, outlier.size=0.2)+
  ggsci::scale_fill_d3()+
  ggsci::scale_color_d3()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

library(ggthemes)
png("/well/combat/users/vkh192/repertoire/out/subsetsource.png",type=c("cairo"), height = 13, width = 13, units = 'in', res = 400)
Heavy %>%
  ggplot(aes(UMAP_1, UMAP_2, color=new_clusters)) +
  geom_point(shape=19, size=0.5, color="black")+
  geom_point(shape=16,size=0.5)+
  scale_colour_tableau("Tableau 20")+
  theme_void()+
  guides(color = guide_legend(override.aes = list(size=3)))+
  facet_wrap(~Source)
dev.off()

#For renyi we need a clone per row vs sequence per column matrix per repertoire
## For CDI the BcellR clonal abundance calculation is essentially the cluster size from RBR pipeline
getCDI <- function(input, B = 10000, conf.level = 0.95) {
  data <- input %>% select(clone_per_replicate, sequence) %>% unique #Here unique B cells will consider all unique BCRs
  counts <- table(data %>% select(clone_per_replicate))
  ID <- input %>% select(COMBAT_ID_Time, Age, Sex ,Source, new_clusters) %>% unique
  options(`mosaic:parallelMessage` = FALSE)
  gini.boot <- mosaic::do(B) * Gini(resample(counts, 20, replace=TRUE))
  names(gini.boot) <- "result"
  boot.mean <- mean(~result, data = gini.boot)
  boot.se <- sd(~result, data = gini.boot)
  boot.median <- median(~result, data = gini.boot)
  g <- as.numeric(1 - conf.level)
  boot.lci <- mosaic::qdata(~result, g/2,data = gini.boot)
  boot.uci <- mosaic::qdata(~result,1-g, data = gini.boot)
  boot.stats <- data.table::data.table(mean=boot.mean, se=boot.se, media=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
  boot.stats <- cbind(ID, boot.stats)
  bootdata <- data.table::data.table(ID, gini_boot=gini.boot)
  return(bootdata)
}

CDI <- mclapply(IDlist, function(x){getCDI(input=x)}, mc.cores=numCores)
CDI <- data.table::rbindlist(CDI)

CDImean <- aggregate(gini_boot.result ~ COMBAT_ID_Time+Source+new_clusters, data=CDI, mean) %>%
  rename(mean_cdi=gini_boot.result)
head(CDImean, 20)

CDImean$Source <- factor(CDImean$Source, levels = c("HV", "Flu","Sepsis", "COVID_LDN","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

png("/well/combat/users/vkh192/repertoire/out/subset_cdi.png",type=c("cairo"), height = 5, width = 8, units = 'in', res = 400)
CDImean %>%
  ggplot(aes(new_clusters, mean_cdi, fill=Source, color=Source))+
  geom_boxplot(alpha=0.4, outlier.size=0.2)+
  ggsci::scale_fill_d3()+
  ggsci::scale_color_d3()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

nrow(CDImean)
CDImean %>% select(COMBAT_ID_Time, Source, new_clusters) %>% duplicated %>% summary
nrow(CEImean)

mergeclonality <- CEImean %>% left_join(CDImean, by=c("COMBAT_ID_Time", "Source", "new_clusters"))
head(mergeclonality, 20)
nrow(mergeclonality)
#%isotype usage plot
png("/well/combat/users/vkh192/repertoire/out/subset_iso.png",type=c("cairo"), height = 5, width = 8, units = 'in', res = 400)
Heavy %>%
  mutate(c_gene= gsub("\\*.*", "", c_gene)) %>%
  filter(c_gene != "None") %>%
  ggplot(aes(new_clusters, fill=c_gene, color=c_gene))+
  geom_bar(position="fill")+
  ggsci::scale_fill_jco()+
  ggsci::scale_color_jco()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  facet_wrap(~Source)
dev.off()

#filter out non-assigned c_genes and create proportion table
Heavyfilt <- Heavy %>% 
  mutate(c_gene= gsub("\\*.*", "", c_gene)) %>%
  filter(c_gene != "None") 

#Frequency table per sample, per cluster, per isotype
Iso_freq<- Heavyfilt %>%
  group_by(COMBAT_ID_Time, Source, new_clusters, c_gene) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

Iso_freq <- spread(Iso_freq %>% select(-n), c_gene, freq)
Iso_freq[is.na(Iso_freq)] <- 0

Iso_freq %>% arrange(COMBAT_ID_Time)
#median SMH
png("/well/combat/users/vkh192/repertoire/out/subset_smh.png",type=c("cairo"), height = 5, width = 8, units = 'in', res = 400)
Heavy %>%
  ggplot(aes(new_clusters,total_mut, fill=Source, color=Source))+
  geom_boxplot(alpha=0.4, outlier.size=0.2)+
  ggsci::scale_fill_d3()+
  ggsci::scale_color_d3()+
  theme(axis.text.x = element_text(angle=90, hjust=1))
dev.off()


#Median SMH per subtype
Heavymut <- aggregate(total_mut ~ COMBAT_ID_Time+Source+new_clusters, data=Heavy, median) %>%
  rename(median_mut_HC=total_mut)
nrow(Heavymut)

#repeat for light chain
Lightmut <- aggregate(total_mut ~ COMBAT_ID_Time+Source+new_clusters, data=Light, median) %>%
  rename(median_mut_LC=total_mut)
head(Lightmut)

merge <- mergeclonality %>% left_join(Iso_freq, by=c("COMBAT_ID_Time", "Source","new_clusters")) %>%
  left_join(Heavymut, by=c("COMBAT_ID_Time", "Source", "new_clusters")) %>% left_join(Lightmut, by=c("COMBAT_ID_Time", "Source", "new_clusters"))
head(merge)
nrow(merge)
#Heavymut <- spread(Heavymut, c_gene, mean_mut_HC)

#clonal sharing
#bootstrap resample clones multiple times, for each bootstrap sampling, perform overlap cdhit clustering. mean the number of overlaps

resampleclones <- function(input, B = 1000, conf.level = 0.95) {
  data <- input[,c("clone_per_replicate")] 
  counts <- (data %>% select(clone_per_replicate))
  ID <- input %>% select(COMBAT_ID_Time, Source, new_clusters) %>% unique
  options(`mosaic:parallelMessage` = FALSE)
  gini.boot <- mosaic::do(B) * (resample(counts, 20, replace=TRUE)) 
  names(gini.boot) <- c("result", "orig.ident", "row", "index")
  bootdata <- data.table::data.table(ID, gini_boot=gini.boot$result)
  return(bootdata)
}

resampclones <- mclapply(IDlist, function(x){resampleclones(input=x)}, mc.cores=numCores)
resampclones <- data.table::rbindlist(resampclones)

#split by combat ID time
cloneslist <- resampclones %>% group_by(COMBAT_ID_Time) %>% group_split()

getoverlap <- function(x){
for(i in unique(x$new_clusters)) {
  ID <- unique(x$COMBAT_ID_Time)
  #get overlaps
  test <- crossprod(table(cloneslist[[1]] %>% ungroup %>% select(gini_boot,new_clusters)))
  test<- (t(test / diag(test))) 
  test <- as.data.frame(test)
  colnames(test) <- paste("vs", colnames(test), sep="_")
  test$new_clusters <- rownames(test)
  test$COMBAT_ID_Time <- ID
}
  return(test)
}

#test <- getoverlap(cloneslist[[1]])

overlap <- mclapply(cloneslist, function(x){getoverlap(x)}, mc.cores=numCores)
overlap <- data.table::rbindlist(overlap)
head(overlap)
#Gather V gene proportions----
#Frequency table per sample, per cluster, per IGHVgene
ighv_freq<- Heavy %>%
  mutate(v_gene_HC= gsub("\\Homsap ", "", v_call)) %>%
  mutate(v_gene_HC= gsub("\\*.*", "", v_gene_HC)) %>%
  mutate(v_gene_HC= gsub("\\-", "_", v_gene_HC))%>%
  group_by(COMBAT_ID_Time, Source, new_clusters, v_gene_HC) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

ighv_freq <- spread(ighv_freq %>% select(-n), v_gene_HC, freq)
ighv_freq[is.na(ighv_freq)] <- 0

#Frequency table per sample, per cluster, per Light V gene
igl_freq<- Light %>%
  mutate(v_gene_LC= gsub("\\Homsap ", "", v_call)) %>%
  mutate(v_gene_LC= gsub("\\*.*", "", v_gene_LC))%>%
  mutate(v_gene_LC= gsub("\\-", "_", v_gene_LC))%>%
  group_by(COMBAT_ID_Time, Source, new_clusters, v_gene_LC) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

igl_freq <- spread(igl_freq %>% select(-n), v_gene_LC, freq)
igl_freq[is.na(igl_freq)] <- 0

#prepare summary metrics for merge
cellcounts<- Heavy %>% 
  select(COMBAT_ID_Time, new_clusters) %>% 
  table %>% 
  reshape2::melt() %>% 
  rename(cells_per_cluster=value)%>% 
  mutate(sampling_depth =ifelse(cells_per_cluster < 5, "Extremely_low", ifelse((cells_per_cluster >=5 & cells_per_cluster< 20), "Low",ifelse((cells_per_cluster >=20 & cells_per_cluster< 100),"Medium" ,ifelse(cells_per_cluster >=100, "High", "Unassigned")))))
nrow(cellcounts)
head(cellcounts)

png("/well/combat/users/vkh192/repertoire/out/cellsperclus.png",type=c("cairo"), height = 5, width = 8, units = 'in', res = 400)
cellcounts %>%
  ggplot(aes(reorder(COMBAT_ID_Time, cells_per_cluster), cells_per_cluster, fill=Source)) +
  geom_col()+
  facet_wrap(~new_clusters)+
  ggsci::scale_fill_d3()+
  theme(axis.text.x=element_blank())
dev.off()

cellcounts %>% group_by(sampling_depth) %>% add_tally() %>% select(sampling_depth, n) %>% unique

#Add Source to cellcounts
Source <- Heavy %>% select(COMBAT_ID_Time, Source) %>% unique
cellcounts <- cellcounts %>% left_join(Source, by=c("COMBAT_ID_Time")) 
#Merge all. To allow for absence of cell populations we will do a full merge
merge2 <- cellcounts %>% 
  full_join(merge%>% select(-Source), by=c("COMBAT_ID_Time","new_clusters")) %>% 
  full_join(overlap, by=c("new_clusters","COMBAT_ID_Time")) %>% 
  full_join(ighv_freq %>% ungroup %>% select(-Source), by=c("new_clusters", "COMBAT_ID_Time")) %>%
  full_join(igl_freq %>% ungroup %>% select(-Source), by=c("new_clusters", "COMBAT_ID_Time"))

head(merge2)
nrow(merge2)

merge2 %>% select(cells_per_cluster) %>% summary
merge2 %>% arrange(COMBAT_ID_Time) %>% head

#Indexing sanity check
heavy_id <- Heavy %>% mutate(ID=paste(COMBAT_ID_Time, Source, sep="_")) %>% select(ID) %>% unique
merge2_id <- merge2 %>% mutate(ID=paste(COMBAT_ID_Time, Source, sep="_")) %>% select(ID) %>% unique 

merge2_id$ID %in% heavy_id$ID %>% summary
heavy_id$ID %in% merge2_id$ID %>% summary

merge2 %>% arrange(COMBAT_ID_Time) %>% head
#Assign removed barcodes as potential doublets from parent database. Run full analysis script until before removal of incomplete VDJs
IMGT$barcode_contig <- paste(IMGT$barcode_seq, IMGT$contig_id, sep="_")
non_full <- IMGT %>% dplyr::filter(v_partial_missing != 0) %>% dplyr::filter(j_partial_missing!=0) %>% select(barcode_contig)

IMGT$parse_chains <- NA

#assign incomplete chains
IMGT <- IMGT %>% mutate(parse_chains=ifelse(!(IMGT$barcode_seq %in% Heavy$barcode_seq), "potential_doublet", "assigned_singlet"))
IMGT <- IMGT %>% mutate(parse_chains=ifelse(barcode_contig %in% non_full$barcode_contig, "incomplete_vdj", parse_chains)) 
head(IMGT)
Parsechainsout <- IMGT %>% select(barcode_seq, COMBAT_ID_Time, locus, parse_chains)
head(Parsechainsout)
Parsechainsout %>% filter(parse_chains=="incomplete_vdj")

write.csv(Parsechainsout, "/well/combat/projects/repertoire/scbcr_table-001/contig_qc-001.csv")
#assign multichain calls
IMGT$barcode_seq %in% Heavy$barcode_seq %>% summary

#Save into shared space
write.csv(merge2, "/well/combat/projects/repertoire/scbcr_table-001/scbcr_table-002.csv")

merge2 %>% filter(new_clusters=="Anergic") %>% arrange(desc(cells_per_cluster)) %>% head
