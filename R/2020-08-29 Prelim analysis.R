#Analysis for COMBAT BCR
setwd("/well/combat/users/vkh192/repertoire/")
.libPaths(c( "~/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0", .libPaths()))
.libPaths(c("/well/combat/users/vkh192/INSTALLED_PROGRAMS/R_MODULES", .libPaths()))
options(repos = 'http://cran.ma.imperial.ac.uk/')

require(tidyverse)
require(patchwork)
require(FSA)
require(ggpubr)
require(gridExtra)
require(grid)
require(cowplot)
require(RColorBrewer)
require(PCAtools)
require(ggsci)
require(rstatix)
require(ggrepel)
require(parallel)
require(UpSetR)
require(Biostrings)

numCores= detectCores()
##Data loading and cleaning ----
clones <- readRDS("/well/combat/users/vkh192/repertoire/data/new_clones.rds")

Heavy <- clones$Heavy
Light <- readRDS("/well/combat/users/vkh192/repertoire/data/Light_chain.rds")

##Clean gene calls by removing alleles and "Homsap" annotations from IMGT, also add new columns with concatenated VJ
Heavy <- Heavy %>% 
  mutate(v_gene_HC= gsub("\\Homsap ", "", v_call)) %>%
  mutate(v_gene_HC= gsub("\\*.*", "", v_gene_HC)) %>%
  mutate(j_gene_HC= gsub("\\Homsap ", "", j_call)) %>%
  mutate(j_gene_HC= gsub("\\*.*", "", j_gene_HC)) %>%
  mutate(VJ_HC= paste(v_gene_HC, j_gene_HC, sep="|")) %>%
  mutate(c_gene= gsub("\\*.*", "", c_gene)) %>%
  mutate(junction_length_HC=nchar(junction_aa)) %>%
  mutate(COMBAT_ID= as.character(COMBAT_ID)) %>%
  mutate(COMBAT_ID_Time= as.character(COMBAT_ID_Time))

Heavy$Source <- factor(Heavy$Source, levels = c("HV", "Flu","Sepsis", "COVID_LDN","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

  
Light <- Light %>% 
  mutate(v_gene_LC= gsub("\\Homsap ", "", v_call)) %>%
  mutate(v_gene_LC= gsub("\\*.*", "", v_gene_LC)) %>%
  mutate(j_gene_LC= gsub("\\Homsap ", "", j_call)) %>%
  mutate(j_gene_LC= gsub("\\*.*", "", j_gene_LC)) %>%
  mutate(VJ_LC= paste(v_gene_LC, j_gene_LC, sep="|")) %>%
  mutate(c_gene= gsub("\\*.*", "", c_gene)) %>%
  mutate(junction_length_LC=nchar(junction_aa))%>%
  mutate(COMBAT_ID= as.character(COMBAT_ID)) %>%
  mutate(COMBAT_ID_Time= as.character(COMBAT_ID_Time))

Light$Source <- factor(Light$Source, levels = c("HV", "Flu","Sepsis", "COVID_LDN","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

##Sample filtering
#Visualise BCR counts and annotations for timepoints and priority
Heavy %>% 
  filter(!grepl("LDN", Source)) %>%
  group_by(COMBAT_ID_Time, Source) %>%
  add_tally() %>%
  select(COMBAT_ID_Time, Source, n, sample_priority, Timepoint) %>%
  unique %>%
  as.data.frame() %>%
  arrange(desc(COMBAT_ID_Time))

#Here I can see that the lowest Timepoint value corresponds to first sampling. 
#For this analysis I will populate metrics for all samples with BCRs >30.
#For comparisons I will compare across first sampling timepoints

print(Heavy %>% 
  group_by(COMBAT_ID_Time, Source) %>%
  add_tally(name="n_bcr") %>%
  filter(n_bcr>30)%>%
  ungroup %>%
  filter(!grepl("LDN", Source)) %>%
  group_by(COMBAT_ID) %>% 
  slice(which.min(Timepoint)) %>%
  ungroup %>%
  group_by(Source) %>%
  add_tally() %>%
  select(Source, n) %>%
  unique %>%
  as.data.frame() %>%
  arrange(desc(Source)))
#After applying proposed filters we have a reasonable distribution of samples. 
#COVID_SEV seems to have the most samples. We will bear this in mind during the analysis. 
#Statistical methods that are robust to sampling depth should be used where possible

#Subset on parameters
samplenames <- Heavy %>% 
  group_by(COMBAT_ID_Time, Source) %>%
  add_tally(name="n_bcr") %>%
  filter(n_bcr>30)%>%
  ungroup() %>%
  filter(!grepl("LDN", Source)) %>%
  group_by(COMBAT_ID) %>% 
  slice(which.min(Timepoint)) %>%
  ungroup() %>%
  select(COMBAT_ID_Time)
#n=117 samples passed filter

Heavy <- Heavy[c(Heavy$COMBAT_ID_Time %in% samplenames$COMBAT_ID_Time), ]
Heavy$Source <- factor(Heavy$Source, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))
Light <- Light[c(Light$COMBAT_ID_Time %in% samplenames$COMBAT_ID_Time), ]
Light$Source <- factor(Light$Source, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))
source <- Heavy %>% ungroup %>% select(COMBAT_ID_Time, Source) %>% unique

##IGHV gene analysis ----
#create proportion table of Heavy V genes. Here each V gene usage is represented as the proportion of BCRs using that V gene in an individuals whole repertoire
HeavyV <- as.data.frame.matrix(table(Heavy$COMBAT_ID_Time, Heavy$v_gene_HC) %>% 
                                 prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.)) %>%
  left_join(source, by=c("COMBAT_ID_Time"))
#Reorder table for kruskal wallis test
HeavyV <- HeavyV %>% select(-COMBAT_ID_Time) %>%
  select(Source, everything())
colnames(HeavyV) <- gsub("-", "_", colnames(HeavyV))

#Kruskal Wallis grouped by Source
ighv_gene_kw<- data.frame(v_gene=NA, p_value=NA, method=NA)[c(1:length(HeavyV)), ]
for (i in 2:ncol(HeavyV))
{
  formula <- as.formula(paste(colnames(HeavyV)[i], " ~ Source", sep=""))
  model <- kruskal.test(formula, data = HeavyV)
  
  ighv_gene_kw$v_gene[i]<- paste(colnames(HeavyV)[i], sep="")
  ighv_gene_kw$p_value[i]<- model$p.value
  ighv_gene_kw$method[i]<- model$method
}
print(ighv_gene_kw %>% arrange(p_value))

ighv_kw_sig <- ighv_gene_kw %>% filter(p_value < 0.05)

#Post hoc Dunn
ighv_gene_dunn <- list()
for (i in 2:ncol(HeavyV))
{
  formula <- as.formula(paste(colnames(HeavyV)[i], " ~ Source", sep=""))
  model <- dunnTest(formula, data = HeavyV, method="bh") 
  result <- model$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = colnames(HeavyV)[i]) %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))
  ighv_gene_dunn[[i]] <- result
}
ighv_gene_dunn<- data.table::rbindlist(ighv_gene_dunn)
ighv_gene_dunn <- ighv_gene_dunn %>% add_x_position(x=".y.", dodge=(0.8))
ighv_dunn_sig <- ighv_gene_dunn %>% filter(P.adj<0.05)
print(ighv_dunn_sig)

p1<- HeavyV %>%
  ungroup() %>%
  reshape2::melt() %>%
  rename(v_gene=variable, group1=Source) %>%
  ggplot(aes(v_gene, value, fill=group1))+
  geom_boxplot(outlier.shape=NA, position = position_dodge(0.8))+
  ggsci::scale_fill_d3()+
  ylim(c(0,0.17))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())+
  ylab("Proportion of repertoire")

png("/well/combat/users/vkh192/repertoire/out/combat_v_gene.png",type=c("cairo"), height = 10, width = 25, units = 'in', res = 400)
p1 / tableGrob(ighv_dunn_sig)
dev.off()

##IGK/LV gene analysis ----
#create proportion table of Heavy V genes. Here each V gene usage is represented as the proportion of BCRs using that V gene in an individuals whole repertoire
LightV <- as.data.frame.matrix(table(Light$COMBAT_ID_Time, Light$v_gene_LC) %>% 
                                 prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.)) %>%
  left_join(source, by=c("COMBAT_ID_Time"))
#Reorder table for kruskal wallis test
LightV <- LightV %>% select(-COMBAT_ID_Time) %>%
  select(Source, everything())
colnames(LightV) <- gsub("-", "_", colnames(LightV))

#Kruskal Wallis grouped by Source
iglv_gene_kw<- data.frame(v_gene=NA, p_value=NA, method=NA)[c(1:length(LightV)), ]
for (i in 2:ncol(LightV))
{
  formula <- as.formula(paste(colnames(LightV)[i], " ~ Source", sep=""))
  model <- kruskal.test(formula, data = LightV)
  
  iglv_gene_kw$v_gene[i]<- paste(colnames(LightV)[i], sep="")
  iglv_gene_kw$p_value[i]<- model$p.value
  iglv_gene_kw$method[i]<- model$method
}
print(iglv_gene_kw %>% arrange(p_value))

iglv_kw_sig <- iglv_gene_kw %>% filter(p_value < 0.05)

#Post hoc Dunn
iglv_gene_dunn <- list()
for (i in 2:ncol(LightV))
{
  formula <- as.formula(paste(colnames(LightV)[i], " ~ Source", sep=""))
  model <- dunnTest(formula, data = LightV, method="bh") 
  result <- model$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = colnames(LightV)[i]) %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))
  iglv_gene_dunn[[i]] <- result
}
iglv_gene_dunn<- data.table::rbindlist(iglv_gene_dunn)
iglv_gene_dunn <- iglv_gene_dunn %>% add_x_position(x=".y.", dodge=(0.8))
iglv_dunn_sig <- iglv_gene_dunn %>% filter(P.adj<0.05)
print(iglv_dunn_sig)

p2<- LightV %>%
  ungroup() %>%
  reshape2::melt() %>%
  rename(v_gene=variable, group1=Source) %>%
  ggplot(aes(v_gene, value, fill=group1))+
  geom_boxplot(outlier.shape=NA, position = position_dodge(0.8))+
  ggsci::scale_fill_d3()+
  ylim(c(0,0.17))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())+
  ylab("Proportion of repertoire")

png("/well/combat/users/vkh192/repertoire/out/combat_iglv_gene.png",type=c("cairo"), height = 15, width = 25, units = 'in', res = 400)
p2 / tableGrob(iglv_dunn_sig)
dev.off()


##IGHV-J pair gene analysis ----
#Resample with replacement uniformly. Subsampling size: n=30, iterations= 100,000
IGHVJ <- Heavy %>% select(COMBAT_ID_Time, VJ_HC)
#IGHVJlist<- split(IGHVJ, IGHVJ$COMBAT_ID_Time)

#Need to write a nicer bootstrap that isnt a loop within a loop!
#IGHVJlist<- lapply(IGHVJlist, function(x){
  #nam <- x[,1] %>% unique()
  #x <- lapply(1:100000, function(i) sample(x$VJ_HC, 30, replace = T))
  #y <- tibble(COMBAT_ID_Time = paste0(nam), VJ_HC=unlist(x))
  #return(y)
  #})
#IGHVJlist[2] %>% head   
#IGHVJ <- data.table::rbindlist(IGHVJlist)

#Create count matrix
HeavyVJ <- as.data.frame.matrix(table(IGHVJ$COMBAT_ID_Time, IGHVJ$VJ_HC) %>%
  prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.))
head(HeavyVJ)

#plot Heavy VJ proportions as boxplot
HeavyVJ <- HeavyVJ %>% 
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  select(-COMBAT_ID_Time) %>%
  select(Source, everything())
colnames(HeavyVJ) <- gsub("-", "_", colnames(HeavyVJ))
colnames(HeavyVJ) <- gsub("\\|", "_", colnames(HeavyVJ))

ncol(HeavyVJ)
#Kruskal Wallis grouped by Source
hvj_gene_kw<- data.frame(v_gene=NA, p_value=NA, method=NA)[c(1:length(HeavyVJ)), ]
for (i in 2:ncol(HeavyVJ))
{
  formula <- as.formula(paste(colnames(HeavyVJ)[i], " ~ Source", sep=""))
  model <- kruskal.test(formula, data = HeavyVJ)
  
  hvj_gene_kw$v_gene[i]<- paste(colnames(HeavyVJ)[i], sep="")
  hvj_gene_kw$p_value[i]<- model$p.value
  hvj_gene_kw$method[i]<- model$method
}
print(hvj_gene_kw %>% arrange(p_value))

hvj_kw_sig <- hvj_gene_kw %>% filter(p_value < 0.05)

#Post hoc Dunn
hvj_gene_dunn <- list()
for (i in 2:ncol(HeavyVJ))
{
  formula <- as.formula(paste(colnames(HeavyVJ)[i], " ~ Source", sep=""))
  model <- dunnTest(formula, data = HeavyVJ, method="bh") 
  result <- model$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = colnames(HeavyVJ)[i]) %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))
   hvj_gene_dunn[[i]] <- result
}
hvj_gene_dunn<- data.table::rbindlist(hvj_gene_dunn)
hvj_gene_dunn <- hvj_gene_dunn %>% add_x_position(x=".y.", dodge=(0.8))

##Extract average log2 fold change. Where one VJ combination is not detectable in another the values are inf. For now I have not corrected this
hvj_gene_dunn$avlog_FC <- NA
for (i in 1:nrow(hvj_gene_dunn)){
  
  hvj_gene_dunn$avlog_FC[i] <- log2(mean(na.omit(HeavyVJ[c(HeavyVJ$Source %in% hvj_gene_dunn$group1[i]),c(colnames(HeavyVJ)%in% hvj_gene_dunn$.y.[i])]))/
                              mean(na.omit(HeavyVJ[c(HeavyVJ$Source %in% hvj_gene_dunn$group2[i]),c(colnames(HeavyVJ)%in% hvj_gene_dunn$.y.[i])])))
}

#select significant values
hvj_dunn_sig <- hvj_gene_dunn %>% filter(P.adj<0.05) 
print(hvj_dunn_sig)
hvj_dunn_sig$group1 <- factor(hvj_dunn_sig$group1, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

png("/well/combat/users/vkh192/repertoire/out/vjpair_dunnmap.png",type=c("cairo"), height = 6, width = 13, units = 'in', res = 400)
hvj_dunn_sig %>%
  ggplot(aes(group2,.y., fill=Z, size=Z, label=round(P.adj, 3))) +
  geom_tile(color="black", size=0.5)+
  geom_text(size=2)+
  facet_wrap(~group1, nrow=1)+
  scale_fill_distiller(palette = 'RdYlBu')+
  #scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  ylab("IGHV-IGHJ gene pair")+
  xlab("")+
  ggtitle("Pairwise Dunn test post Kruskal-Wallis test with BH correction")
dev.off()


png("/well/combat/users/vkh192/repertoire/out/vjpair_dunnz.png",type=c("cairo"), height = 7, width = 6, units = 'in', res = 400)
hvj_dunn_sig %>%
  ggplot(aes(.y., Z, fill=group1, size=P.adj, label=round(P.adj, 3))) +
  geom_point(shape=21, color="gray", alpha=0.7)+
  ggsci::scale_fill_d3()+
  theme_bw()+
  #scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  ylab("Z statistic")+
  xlab("")+
  coord_flip()+
  ggtitle("IGHVJ pairs vs Z statistic of group1 vs.")
dev.off()


##Light VJ analysis ----
#Create prob matrix
LightVJ <- as.data.frame.matrix(table(Light$COMBAT_ID_Time, Light$VJ_LC)%>% 
                                  prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.))

#clean Light VJ annotations to be compatible with code
LightVJ <- LightVJ %>% 
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  select(-COMBAT_ID_Time) %>%
  select(Source, everything())
colnames(LightVJ) <- gsub("-", "_", colnames(LightVJ))
colnames(LightVJ) <- gsub("\\|", "_", colnames(LightVJ))

#Kruskal Wallis grouped by Source
lvj_gene_kw<- data.frame(v_gene=NA, p_value=NA, method=NA)[c(1:length(LightVJ)), ]
for (i in 2:ncol(LightVJ))
{
  formula <- as.formula(paste(colnames(LightVJ)[i], " ~ Source", sep=""))
  model <- kruskal.test(formula, data = LightVJ)
  
  lvj_gene_kw$v_gene[i]<- paste(colnames(LightVJ)[i], sep="")
  lvj_gene_kw$p_value[i]<- model$p.value
  lvj_gene_kw$method[i]<- model$method
}
print(lvj_gene_kw %>% arrange(p_value))

lvj_kw_sig <- lvj_gene_kw %>% filter(p_value < 0.05)

#Post hoc Dunn
lvj_gene_dunn <- list()
for (i in 2:ncol(LightVJ))
{
  formula <- as.formula(paste(colnames(LightVJ)[i], " ~ Source", sep=""))
  model <- dunnTest(formula, data = LightVJ, method="bh") 
  result <- model$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = colnames(LightVJ)[i]) %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))
  lvj_gene_dunn[[i]] <- result
}
lvj_gene_dunn<- data.table::rbindlist(lvj_gene_dunn)
lvj_gene_dunn <- lvj_gene_dunn %>% add_x_position(x=".y.", dodge=(0.8))

##Extract average log2 fold change. Where one VJ combination is not detectable in another the values are inf. For now I have not corrected this
lvj_gene_dunn$avlog_FC <- NA
for (i in 1:nrow(lvj_gene_dunn)){
  
  lvj_gene_dunn$avlog_FC[i] <- log2(mean(na.omit(LightVJ[c(LightVJ$Source %in% lvj_gene_dunn$group1[i]),c(colnames(LightVJ)%in% lvj_gene_dunn$.y.[i])]))/
                                      mean(na.omit(LightVJ[c(LightVJ$Source %in% lvj_gene_dunn$group2[i]),c(colnames(LightVJ)%in% lvj_gene_dunn$.y.[i])])))
}

#select significant values
lvj_dunn_sig <- lvj_gene_dunn %>% filter(P.adj<0.05) 
print(lvj_dunn_sig)
lvj_dunn_sig$group1 <- factor(lvj_dunn_sig$group1, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

png("/well/combat/users/vkh192/repertoire/out/lvjpair_dunnmap.png",type=c("cairo"), height = 6, width = 13, units = 'in', res = 400)
lvj_dunn_sig %>%
  ggplot(aes(group2,.y., fill=Z, size=Z, label=round(P.adj, 3))) +
  geom_tile(color="black", size=0.5)+
  geom_text(size=2)+
  facet_wrap(~group1, nrow=1)+
  scale_fill_distiller(palette = 'RdYlBu')+
  #scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  ylab("Light chain VJ gene pair")+
  xlab("")+
  ggtitle("Pairwise Dunn test post Kruskal-Wallis test with BH correction")
dev.off()


png("/well/combat/users/vkh192/repertoire/out/lvjpair_dunnz.png",type=c("cairo"), height = 7, width = 6, units = 'in', res = 400)
lvj_dunn_sig %>%
  ggplot(aes(.y., Z, fill=group1, size=P.adj, label=round(P.adj, 3))) +
  geom_point(shape=21, color="gray", alpha=0.7)+
  ggsci::scale_fill_d3()+
  theme_bw()+
  #scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  ylab("Z statistic")+
  xlab("")+
  coord_flip()+
  ggtitle("Light chain VJ pairs vs Z statistic of group1 vs.")
dev.off()


##Heavy chain isotype usage and mutations----
#filter out non-assigned c_genes and create proportion table
Heavyfilt <- Heavy %>% 
  filter(c_gene != "None") 
HeavyC <- as.data.frame.matrix(table(Heavyfilt$COMBAT_ID_Time, Heavyfilt$c_gene)%>% 
                                 prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.))

HeavyC <- HeavyC %>% 
  left_join(source, by=c("COMBAT_ID_Time"))

#KW+Dunn
HeavyC %>%
  group_by(Source) %>%
  summarise(N = n()) 
HeavyC <- HeavyC %>% select(-COMBAT_ID_Time) %>%
  select(Source, everything())

c_gene_kw<- data.frame(c_gene=NA, p_value=NA, method=NA)[c(1:length(HeavyC)), ]
for (i in 2:ncol(HeavyC))
{
  formula <- as.formula(paste(colnames(HeavyC)[i], " ~ Source", sep=""))
  model <- kruskal.test(formula, data = HeavyC)
  
  c_gene_kw$c_gene[i]<- paste(colnames(HeavyC)[i], sep="")
  c_gene_kw$p_value[i]<- model$p.value
  c_gene_kw$method[i]<- model$method
}

#Post hoc Dunn
c_gene_dunn <- list()
for (i in 2:ncol(HeavyC))
{
  formula <- as.formula(paste(colnames(HeavyC)[i], " ~ Source", sep=""))
  model <- dunnTest(formula, data = HeavyC, method="bh") 
  result <- model$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = colnames(HeavyC)[i]) %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))
  c_gene_dunn[[i]] <- result
}
c_gene_dunn<- data.table::rbindlist(c_gene_dunn)
c_gene_dunn <- c_gene_dunn %>% mutate(c_gene= .y.) 
#manually add xvalues based on group1 and 2
c_gene_dunn <- c_gene_dunn %>% 
  mutate(xmin = ifelse(grepl("HV", group1), 1, ifelse(grepl("Flu", group1), 2, ifelse(grepl("Sepsis", group1), 3, ifelse(grepl("HCW", group1),
  4, ifelse(grepl("COVID_MILD", group1), 5, ifelse(grepl("SEV", group1), 6, 7))))))) %>%
  mutate(xmax = ifelse(grepl("HV", group2), 1, ifelse(grepl("Flu", group2), 2, ifelse(grepl("Sepsis", group2), 3, ifelse(grepl("HCW", group2),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group2), 5, ifelse(grepl("SEV", group2), 6, 7)))))))

##Extract average log2 fold change. Where one VJ combination is not detectable in another the values are inf. For now I have not corrected this
c_gene_dunn$avlog_FC <- NA
for (i in 1:nrow(c_gene_dunn)){
  
  c_gene_dunn$avlog_FC[i] <- log2(mean(na.omit(HeavyC[c(HeavyC$Source %in% c_gene_dunn$group1[i]),c(colnames(HeavyC)%in% c_gene_dunn$.y.[i])]))/
                                      mean(na.omit(HeavyC[c(HeavyC$Source %in% c_gene_dunn$group2[i]),c(colnames(HeavyC)%in% c_gene_dunn$.y.[i])])))
}

#select significant values
c_dunn_sig <- c_gene_dunn %>% filter(P.adj<0.05) 
print(c_dunn_sig)
c_dunn_sig$group1 <- factor(c_dunn_sig$group1, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

c_gene_kw %>% filter(p_value < 0.05)

png("/well/combat/users/vkh192/repertoire/out/combat_iso.png",type=c("cairo"), height = 5, width = 9, units = 'in', res = 400)
HeavyC %>% 
  reshape2::melt(id=c("Source")) %>%
  rename(c_gene=variable, group1=Source) %>%
  ggplot(aes(group1, value, fill=group1))+
  geom_violin(alpha=0.4)+
  geom_jitter(width=0.2, shape=21, alpha=0.4, size=0.8)+
  stat_summary(fill="white", fun.y=mean, geom="point", shape=23, size=0.8)+
  facet_wrap(~c_gene, scales = "free_y")+
  stat_pvalue_manual(c_gene_dunn, label = "p.signif",
                     y.position = c(0.6, 0.66, 0.72, 0.25, 0.30, 0.35, 0.40, 0.0025, 0.0030, 0.0035, 0.004,0.85, 0.90, 0.95, 1, 1.05), hide.ns = TRUE)+
  ylab("Proportion of repertoire")+
  ggsci::scale_fill_d3()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  xlab("")
dev.off()


#Heavy iso mutations
#Create median total mutations table per individual repertoire. We choose median here as mutations often follow a poisson distribution.
Heavymut <- aggregate(total_mut ~ c_gene+COMBAT_ID_Time+Source, data=Heavyfilt, median) %>%
  rename(mean_mut_HC=total_mut)

Heavymut <- spread(Heavymut, c_gene, mean_mut_HC)
Heavymut[is.na(Heavymut)] <- 0

#Remove IGHE as there is only one observation and it causes the KW loop to fail
Heavymut <- Heavymut %>% select(-IGHE)

mut_kw<- data.frame(c_gene=NA, p_value=NA, method=NA)[c(1:length(Heavymut)), ]
for (i in 3:ncol(Heavymut))
{
  formula <- as.formula(paste(colnames(Heavymut)[i], " ~ Source", sep=""))
  model <- kruskal.test(formula, data = as.matrix(Heavymut))
  mut_kw$c_gene[i]<- paste(colnames(Heavymut)[i], sep="")
  mut_kw$p_value[i]<- model$p.value
  mut_kw$method[i]<- model$method
}

#Post hoc Dunn
mut_dunn <- list()
for (i in 3:ncol(Heavymut))
{
  formula <- as.formula(paste(colnames(Heavymut)[i], " ~ Source", sep=""))
  model <- dunnTest(formula, data = Heavymut, method="bh") 
  result <- model$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = colnames(Heavymut)[i]) %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))
  mut_dunn[[i]] <- result
}
mut_dunn<- data.table::rbindlist(mut_dunn)
mut_dunn <- mut_dunn %>% mutate(c_gene= .y.) 
#manually add xvalues based on group1 and 2
mut_dunn <- mut_dunn %>% 
  mutate(xmin = ifelse(grepl("HV", group1), 1, ifelse(grepl("Flu", group1), 2, ifelse(grepl("Sepsis", group1), 3, ifelse(grepl("HCW", group1),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group1), 5, ifelse(grepl("SEV", group1), 6, 7))))))) %>%
  mutate(xmax = ifelse(grepl("HV", group2), 1, ifelse(grepl("Flu", group2), 2, ifelse(grepl("Sepsis", group2), 3, ifelse(grepl("HCW", group2),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group2), 5, ifelse(grepl("SEV", group2), 6, 7)))))))

##Extract average log2 fold change. Where one VJ combination is not detectable in another the values are inf. For now I have not corrected this
mut_dunn$avlog_FC <- NA
for (i in 1:nrow(mut_dunn)){
  
  mut_dunn$avlog_FC[i] <- log2(mean(na.omit(Heavymut[c(Heavymut$Source %in% mut_dunn$group1[i]),c(colnames(Heavymut)%in% mut_dunn$.y.[i])]))/
                                    mean(na.omit(Heavymut[c(Heavymut$Source %in% mut_dunn$group2[i]),c(colnames(Heavymut)%in% mut_dunn$.y.[i])])))
}

#select significant values
mut_dunn_sig <- mut_dunn %>% filter(P.adj<0.05) 
print(mut_dunn_sig)
mut_dunn$group1 <- factor(mut_dunn$group1, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

mut_kw %>% filter(p_value < 0.05)



head(Heavy$Source)
png("/well/combat/users/vkh192/repertoire/out/combat_mut.png",type=c("cairo"), height = 5, width = 9, units = 'in', res = 400)
Heavymut %>% 
  reshape2::melt() %>%
  rename(c_gene=variable, group1=Source) %>%
  filter(c_gene != "None") %>%
  ggplot(aes(group1, value, fill=group1))+
  geom_violin(alpha=0.5)+
  geom_jitter(width=0.2, shape=21, alpha=0.4, size=0.8)+
  stat_summary(fill="white", fun.y=mean, geom="point", shape=23, size=0.8)+
  stat_pvalue_manual(mut_dunn, label = "p.signif",
                     y.position = c(70,74,78,82,86,90,94,50,30,30), hide.ns = TRUE)+  
  ggsci::scale_fill_d3()+
  facet_wrap(~c_gene)+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  ylab("Median total mutations")+
  xlab("")
dev.off()


##Light KL use **need to do---- 
LightC <- as.data.frame.matrix(table(Light$COMBAT_ID_Time, Light$c_gene)%>% 
                                 prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.))





##CDR3 lengths ----

Heavycdr3len <- aggregate(junction_length_HC ~ c_gene+COMBAT_ID_Time+Source, data=Heavyfilt, mean) %>%
  rename(mean_junclen_HC=junction_length_HC)

Heavycdr3len <- spread(Heavycdr3len, c_gene, mean_junclen_HC)


#Remove IGHE as there is only one observation and it causes the KW loop to fail
Heavycdr3len <- Heavycdr3len %>% select(-IGHE)

cdr3_kw<- data.frame(c_gene=NA, p_value=NA, method=NA)[c(1:length(Heavycdr3len)), ]
for (i in 3:ncol(Heavycdr3len))
{
  formula <- as.formula(paste(colnames(Heavycdr3len)[i], " ~ Source", sep=""))
  model <- kruskal.test(formula, data = as.matrix(Heavycdr3len))
  cdr3_kw$c_gene[i]<- paste(colnames(Heavycdr3len)[i], sep="")
  cdr3_kw$p_value[i]<- model$p.value
  cdr3_kw$method[i]<- model$method
}

#Post hoc Dunn
cdr_dunn <- list()
for (i in 3:ncol(Heavycdr3len))
{
  formula <- as.formula(paste(colnames(Heavycdr3len)[i], " ~ Source", sep=""))
  model <- dunnTest(formula, data = Heavycdr3len, method="bh") 
  result <- model$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = colnames(Heavycdr3len)[i]) %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))
  cdr_dunn[[i]] <- result
}
cdr_dunn<- data.table::rbindlist(cdr_dunn)
cdr_dunn <- cdr_dunn %>% mutate(c_gene= .y.) 
#manually add xvalues based on group1 and 2
cdr_dunn <- cdr_dunn %>% 
  mutate(xmin = ifelse(grepl("HV", group1), 1, ifelse(grepl("Flu", group1), 2, ifelse(grepl("Sepsis", group1), 3, ifelse(grepl("HCW", group1),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group1), 5, ifelse(grepl("SEV", group1), 6, 7))))))) %>%
  mutate(xmax = ifelse(grepl("HV", group2), 1, ifelse(grepl("Flu", group2), 2, ifelse(grepl("Sepsis", group2), 3, ifelse(grepl("HCW", group2),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group2), 5, ifelse(grepl("SEV", group2), 6, 7)))))))

##Extract average log2 fold change. Where one VJ combination is not detectable in another the values are inf. For now I have not corrected this
cdr_dunn$avlog_FC <- NA
for (i in 1:nrow(mut_dunn)){
  
  cdr_dunn$avlog_FC[i] <- log2(mean(na.omit(Heavycdr3len[c(Heavycdr3len$Source %in% cdr_dunn$group1[i]),c(colnames(Heavycdr3len)%in% cdr_dunn$.y.[i])]))/
                                 mean(na.omit(Heavycdr3len[c(Heavycdr3len$Source %in% cdr_dunn$group2[i]),c(colnames(Heavycdr3len)%in% cdr_dunn$.y.[i])])))
}

#select significant values
cdr_dunn_sig <- cdr_dunn %>% filter(P.adj<0.05) 
print(cdr_dunn_sig)
cdr_dunn$group1 <- factor(cdr_dunn$group1, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

cdr_kw %>% filter(p_value < 0.05)


#plot mean cdr3 lengths
png("/well/combat/users/vkh192/repertoire/out/combat_cdr3.png",type=c("cairo"), height = 5, width = 9, units = 'in', res = 400)
Heavycdr3len %>%
  reshape2::melt() %>%
  rename(c_gene=variable, group1=Source) %>%
  filter(c_gene != "None") %>%
  ggplot(aes(x = group1, y = value, fill=group1)) +
  geom_violin(alpha=0.5)+
  geom_jitter(width=0.2, shape=21, alpha=0.4, size=0.8)+
  stat_summary(fill="white", fun.y=mean, geom="point", shape=23, size=0.8)+
    stat_pvalue_manual(cdr_dunn, label = "p.signif",
                       y.position = c(25), hide.ns = TRUE)+  
    ggsci::scale_fill_d3()+
    facet_wrap(~c_gene)+
  ylab("Mean cdr3 length")+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()


##junction hydrophobicity need to do----

Heavyhydro <- aggregate(junction_hydrophobicity ~ COMBAT_ID_Time, data=Heavy, median) %>%
  rename(hydro_HC=junction_hydrophobicity)

Lighthydro <- aggregate(junction_hydrophobicity ~ COMBAT_ID_Time, data=Light, mean) %>%
  rename(hydro_junclen_LC=junction_hydrophobicity)

##Get clonality measures----

#Clonal expansion index is the simpsons index for each unique VDJ identity with uniform subsampling to the same number of unique BCRs

IDlist<- Heavy %>% 
  group_by(COMBAT_ID_Time) %>%
  group_split()

#number of unique vdjs by nucleotide
print(ID1 %>% select(sequence) %>% unique %>% nrow)

library(vegan)
library(mosaic)
library(DescTools)

#For CEI, we must first group by barcodes and VDJ sequence per sample to collapse down to unique. Resampling will be done to the lowest sample size
Heavy %>% 
  select(COMBAT_ID_Time, sequence) %>%
  unique %>% 
  group_by(COMBAT_ID_Time) %>% 
  add_tally() %>%
  ungroup%>%
  select(COMBAT_ID_Time, n)%>% 
  unique %>%
  as.data.frame() %>%
  arrange(n)
#Here we see the smallest number of individual BCRs is 25, therefore subsampling must capture at least this many. 
#If we replace during subsampling, we can still estimate diversity indirectly vs. sampling without replacement which only captures the same diversity

getCEI <- function(input, B = 1000 , conf.level = 0.95) {
  data <- input[,c("sequence")]
  counts <- table(data)
  ID <- input %>% select(COMBAT_ID_Time, Age, Sex ,Source) %>% unique
  options(`mosaic:parallelMessage` = FALSE)
  gini.boot <- mosaic::do(B) * Gini(resample(counts, 25, replace=TRUE))
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

#For renyi we need a clone per row vs sequence per column matrix per repertoire
## For CDI the BcellR clonal abundance calculation is essentially the cluster size from RBR pipeline
getCDI <- function(input, B = 1000, conf.level = 0.95) {
  data <- input %>% select(clone_per_replicate, sequence) %>% unique #Here unique B cells will consider all unique BCRs
  counts <- table(data %>% select(clone_per_replicate))
  ID <- input %>% select(COMBAT_ID_Time, Age, Sex ,Source) %>% unique
  options(`mosaic:parallelMessage` = FALSE)
  gini.boot <- mosaic::do(B) * Gini(resample(counts, 25, replace=TRUE))
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

#Statistics for CEI----
CEImean <- aggregate(gini_boot.result ~ COMBAT_ID_Time+Source, data=CEI, mean) %>%
  rename(mean_cei=gini_boot.result)
CEImean

#KW
cei_kw <- kruskal.test(mean_cei~Source, data = (CEImean))


#Post hoc Dunn
cei_dunn <- dunnTest(mean_cei~Source, data = CEImean, method="bh") 
cei_dunn <- cei_dunn$res %>% 
    as.data.frame()%>%
    separate(Comparison, c("group1", "group2"), "\\ - ") %>%
    mutate(.y. = "mean_cei") %>%
    mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                            ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                   ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                          ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                                 ifelse(P.adj <=0.0001, "****", "na"))))))


#manually add xvalues based on group1 and 2
cei_dunn <- cei_dunn %>% 
  mutate(xmin = ifelse(grepl("HV", group1), 1, ifelse(grepl("Flu", group1), 2, ifelse(grepl("Sepsis", group1), 3, ifelse(grepl("HCW", group1),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group1), 5, ifelse(grepl("SEV", group1), 6, 7))))))) %>%
  mutate(xmax = ifelse(grepl("HV", group2), 1, ifelse(grepl("Flu", group2), 2, ifelse(grepl("Sepsis", group2), 3, ifelse(grepl("HCW", group2),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group2), 5, ifelse(grepl("SEV", group2), 6, 7)))))))

#select significant values
cei_dunn_sig <- cei_dunn %>% filter(P.adj<0.05) 
print(cei_dunn_sig)
cei_dunn$group1 <- factor(cei_dunn$group1, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

#plot with sig values
png("/well/combat/users/vkh192/repertoire/out/combat_cei.png",type=c("cairo"), height = 4, width = 5, units = 'in', res = 400)
CEImean %>%
  rename(group1=Source) %>%
  ggplot(aes(group1, mean_cei, fill=group1))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  geom_jitter(width=0.2, shape=21, alpha=0.4, size=0.8)+
  stat_pvalue_manual(cei_dunn, label = "p.signif",
                     y.position = c(0.37, 0.41, 0.45, 0.49), hide.ns = TRUE)+  
  ggsci::scale_fill_d3()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()





#Statistics for CDI----
CDImean <- aggregate(gini_boot.result ~ COMBAT_ID_Time+Source, data=CDI, mean) %>%
  rename(mean_cdi=gini_boot.result)
CDImean

#KW
cdi_kw <- kruskal.test(mean_cdi~Source, data = (CDImean))


#Post hoc Dunn
cdi_dunn <- dunnTest(mean_cdi~Source, data = CDImean, method="bh") 
cdi_dunn <- cdi_dunn$res %>% 
  as.data.frame()%>%
  separate(Comparison, c("group1", "group2"), "\\ - ") %>%
  mutate(.y. = "mean_cei") %>%
  mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                          ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                 ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                        ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                               ifelse(P.adj <=0.0001, "****", "na"))))))


#manually add xvalues based on group1 and 2
cdi_dunn <- cdi_dunn %>% 
  mutate(xmin = ifelse(grepl("HV", group1), 1, ifelse(grepl("Flu", group1), 2, ifelse(grepl("Sepsis", group1), 3, ifelse(grepl("HCW", group1),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group1), 5, ifelse(grepl("SEV", group1), 6, 7))))))) %>%
  mutate(xmax = ifelse(grepl("HV", group2), 1, ifelse(grepl("Flu", group2), 2, ifelse(grepl("Sepsis", group2), 3, ifelse(grepl("HCW", group2),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group2), 5, ifelse(grepl("SEV", group2), 6, 7)))))))

#select significant values
cdi_dunn_sig <- cdi_dunn %>% filter(P.adj<0.05) 
print(cdi_dunn_sig)
cdi_dunn$group1 <- factor(cdi_dunn$group1, levels = c("HV", "Flu","Sepsis","COVID_HCW_MILD", "COVID_MILD", "COVID_SEV", "COVID_CRIT"))

png("/well/combat/users/vkh192/repertoire/out/combat_cdi.png",type=c("cairo"), height = 4, width = 5, units = 'in', res = 400)
CDImean %>%
  rename(group1=Source) %>%
  ggplot(aes(group1, mean_cdi, fill=group1))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  geom_jitter(width=0.2, shape=21, alpha=0.4, size=0.8)+
  stat_pvalue_manual(cdi_dunn, label = "p.signif",
                     y.position = c(0.32, 0.36, 0.40, 0.44, 0.48, 0.52), hide.ns = TRUE)+  
  ggsci::scale_fill_d3()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

##Shannons diversity ----
#Shannon function
shannon <- function(data) {
  counts <- table(data)
  p.i <- counts/sum(counts)
  H <- (-1) * sum(p.i * log(p.i))
  return(H)
}

#Bootstrap shannon
boot.shannon <- function(input, B = 1000, conf.level = 0.95) {
  data <- input$clone_per_replicate
  ID <- input %>% select(COMBAT_ID_Time, Age, Sex ,Source) %>% unique
  shannon.boot <- mosaic::do(B) * shannon(resample(data, 30, replace=TRUE))
  names(shannon.boot) <- "result"
  boot.mean <- mean(~result, data = shannon.boot)
  boot.se <- sd(~result, data = shannon.boot)
  boot.median <- median(~result, data = shannon.boot)
  g <- as.numeric(1 - conf.level)
  boot.lci <- mosaic::qdata(~result, g/2,data = shannon.boot)
  boot.uci <- mosaic::qdata(~result,1-g, data = shannon.boot)
  boot.stats <- data.table::data.table(mean=boot.mean, se=boot.se, media=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
  boot.stats <- cbind(ID, boot.stats)
  bootdata <- data.table::data.table(ID, shannon_boot=shannon.boot)
  return(bootdata)
}

#Get shannon's diversity
shannond <- mclapply(IDlist, function(x){boot.shannon(x)}, mc.cores=numCores)
shannond <- data.table::rbindlist(shannond)

#Statistics for shannons diversity----
shanmean <- aggregate(shannon_boot.result ~ COMBAT_ID_Time+Source, data=shannond, mean) %>%
  rename(mean_shan=shannon_boot.result)
shanmean

#Remove IGHE as there is only one observation and it causes the KW loop to fail
shan_kw <- kruskal.test(mean_shan~Source, data = (shanmean))


#Post hoc Dunn
shan_dunn <- dunnTest(mean_shan~Source, data = shanmean, method="bh") 
shan_dunn <- shan_dunn$res %>% 
  as.data.frame()%>%
  separate(Comparison, c("group1", "group2"), "\\ - ") %>%
  mutate(.y. = "mean_cei") %>%
  mutate(p.signif= ifelse(P.adj >=0.05, "ns", 
                          ifelse(P.adj >=0.01 & P.adj <=0.05, "*",
                                 ifelse(P.adj >= 0.001 & P.adj <=0.01, "**",
                                        ifelse(P.adj >= 0.0001 & P.adj <=0.001 , "***",
                                               ifelse(P.adj <=0.0001, "****", "na"))))))


#manually add xvalues based on group1 and 2
shan_dunn <- shan_dunn %>% 
  mutate(xmin = ifelse(grepl("HV", group1), 1, ifelse(grepl("Flu", group1), 2, ifelse(grepl("Sepsis", group1), 3, ifelse(grepl("HCW", group1),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group1), 5, ifelse(grepl("SEV", group1), 6, 7))))))) %>%
  mutate(xmax = ifelse(grepl("HV", group2), 1, ifelse(grepl("Flu", group2), 2, ifelse(grepl("Sepsis", group2), 3, ifelse(grepl("HCW", group2),
                                                                                                                         4, ifelse(grepl("COVID_MILD", group2), 5, ifelse(grepl("SEV", group2), 6, 7)))))))

#select significant values
shan_dunn_sig <- shan_dunn %>% filter(P.adj<0.05) 
print(shan_dunn_sig)

png("/well/combat/users/vkh192/repertoire/out/combat_shan.png",type=c("cairo"), height = 4, width = 5, units = 'in', res = 400)
shanmean %>%
  rename(group1=Source) %>%
  ggplot(aes(group1, mean_shan, fill=group1))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  geom_jitter(width=0.2, shape=21, alpha=0.4, size=0.8)+
  stat_pvalue_manual(shan_dunn, label = "p.signif",
                     y.position = c(3.7, 3.9, 4.1, 4.3, 4.5, 4.7), hide.ns = TRUE)+  
  ggsci::scale_fill_d3()+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()

##Clonal overlap----
#During BcellR clonal assigment step, we undertook a glonal clone assignment. As edgelists are forms originally on the combined dataset, we have a glonal clone ID
#This captures shared clones by similarity and clonal relatedness not just 100% identity
shared <- Heavy %>% select(clone_global, Source) %>% unique %>% group_by(clone_global) %>% add_tally() %>% filter(n>=2) %>% arrange(clone_global)
shared %>%  ungroup %>% select(Source) %>% unique

Heavy[c(Heavy$clone_global %in% shared$clone_global), ] %>% select(Source) %>% unique

VENN.LIST <- Heavy %>% select(Source, clone_global) %>% split(., .$Source, drop = TRUE)
VENN.LIST <- lapply(VENN.LIST, function(x){x <- x$clone_global})

names(VENN.LIST)
png("/well/combat/users/vkh192/repertoire/out/combat_venn.png",type=c("cairo"), height = 5, width = 10, units = 'in', res = 400)
upset(fromList(VENN.LIST), order.by="degree", point.size = 2, line.size = 1, nsets=100)
dev.off()

#Extract most abundant intersects
#COVID MILD vs SEV
shared_cms <- Heavy %>% 
  filter(grepl("COVID_MILD|COVID_SEV", Source)) %>% 
  select(clone_global, Source) %>% 
  unique() %>%
  group_by(clone_global) %>% 
  add_tally() %>% 
  filter(n>=2) %>% 
  arrange(clone_global)

shared_cms %>%  
  ungroup() %>% 
  select(clone_global) %>% 
  unique() %>% 
  as.data.frame

conv_clono<- Heavy[c(Heavy$clone_global %in% shared_cms$clone_global), ] %>% 
  select(gplex,umis, clone_global, Source, v_gene_HC, j_gene_HC, junction_aa, total_mut, c_gene) %>% 
  unique() %>% 
  arrange(clone_global) %>%
  as.data.frame()


##Antigen specific ----
#upload sarscov2 database to rescomp
#prepare OPIG db
covdb <- readr::read_csv("/well/combat/users/vkh192/repertoire/data/CoV-AbDab_270820.csv")

covdb <- covdb %>% filter(grepl("Human", Origin))
covdb$CDRH3 <- paste("C",covdb$CDRH3, "W", sep="")
covdbHC <- covdb %>% mutate(v_gene_HC= gsub("\\ .*", "", `Heavy V Gene`)) %>% 
  mutate(j_gene_HC= gsub("\\ .*", "", `Heavy J Gene`)) %>% 
  dplyr::rename(junction_aa = CDRH3, epitope=`Protein + Epitope`, neutralising = `Neutralising Vs`) %>% 
  mutate(ID = Name) %>%
  mutate(VJ_HC = paste(.$v_gene_HC, .$j_gene_HC, sep="|")) %>%
  dplyr::select(ID, v_gene_HC, j_gene_HC, junction_aa, VJ_HC, epitope, neutralising, Sources)

covdbLC <- covdb %>% mutate(v_gene_LC= gsub("\\ .*", "", `Light V Gene`)) %>% 
  mutate(j_gene_LC= gsub("\\ .*", "", `Light J Gene`)) %>% 
  dplyr::rename(junction_aa_LC = CDRL3) %>% 
  mutate(ID = Name) %>%
  mutate(VJ_LC = paste(v_gene_LC, j_gene_LC, sep="|")) %>%
  select(ID, v_gene_LC,j_gene_LC, junction_aa_LC, VJ_LC) %>%
  mutate(junction_aa_LC = paste("C", junction_aa_LC, "F", sep=""))

#prepare COMBAT seqs
combatvj <- Heavy %>% 
  select(barcode_seq, v_gene_HC, j_gene_HC, junction_aa, VJ_HC, Source, COMBAT_ID, total_mut, clone_per_replicate, clone_global, clonal_abundance, c_gene, name) %>%
  dplyr::rename(ID=barcode_seq)


#load screaton mAbs
heavy_s <- readRDS("/well/combat/users/vkh192/repertoire/data/heavy_s.rds")
light_s <- readRDS("/well/combat/users/vkh192/repertoire/data/light_s.rds")

heavy_s <- heavy_s %>% 
  mutate(v_gene_HC= gsub("\\Homsap ", "", v_call)) %>%
  mutate(v_gene_HC= gsub("\\*.*", "", v_gene_HC)) %>%
  mutate(j_gene_HC= gsub("\\Homsap ", "", j_call)) %>%
  mutate(j_gene_HC= gsub("\\*.*", "", j_gene_HC)) %>%
  mutate(VJ_HC= paste(v_gene_HC, j_gene_HC, sep="|")) %>%
  mutate(junction_length_HC=nchar(junction_aa)) %>%
  dplyr::rename(ID=sample, COMBAT_ID=contig, epitope=barcode) %>%
  select(ID, v_gene_HC, j_gene_HC, junction_aa, VJ_HC, epitope)

light_s <- light_s$Light %>% 
  mutate(v_gene_LC= gsub("\\Homsap ", "", v_call)) %>%
  mutate(v_gene_LC= gsub("\\*.*", "", v_gene_LC)) %>%
  mutate(j_gene_LC= gsub("\\Homsap ", "", j_call)) %>%
  mutate(j_gene_LC= gsub("\\*.*", "", j_gene_LC)) %>%
  mutate(VJ_LC= paste(v_gene_LC, j_gene_LC, sep="|")) %>%
  mutate(junction_length_HC= nchar(junction_aa)) %>%
  dplyr::rename(ID=sample, COMBAT_ID=contig, epitope=barcode, junction_aa_LC=junction_aa) %>%
  select(ID, v_gene_LC, j_gene_LC, junction_aa_LC, VJ_LC)

#Rbind all specific heavy chains + ight chains
head(covdbHC)
head(heavy_s)
heavy_s$neutralising <- NA
heavy_s$Sources <- "Screaton"

specific_HC <- rbind(covdbHC, heavy_s)

head(covdbLC)
head(light_s)
specific_LC <- rbind(covdbLC, light_s)

#cluster by cdr3 levdist

library(stringdist)
spec_cdrdist <- stringdistmatrix(combatvj$junction_aa, specific_HC$junction_aa, method="lv")

rownames(spec_cdrdist) <- as.character(combatvj$ID)
colnames(spec_cdrdist) <- as.character(specific_HC$ID)

spec_cdrdist[upper.tri(spec_cdrdist)] <- NA
spec_cdrdist <- reshape2::melt(spec_cdrdist)
colnames(spec_cdrdist) <- c("ID", "seq_to", "lv_dist")

png("/well/combat/users/vkh192/repertoire/out/cov_disthist.png",type=c("cairo"), height = 6, width = 7, units = 'in', res = 400)
ggplot(spec_cdrdist, aes(lv_dist))+
  geom_histogram(binwidth=1)
dev.off()

#pass threshold of lv=5
match <- spec_cdrdist %>% filter(lv_dist <=2) %>% arrange(desc(lv_dist)) %>% unique

combatvj$ID <- as.character(combatvj$ID)
match$ID <- as.character(match$ID)
match <- match %>% left_join(combatvj, by=c("ID"))
match <- match %>% left_join(specific_HC, by=c("seq_to"="ID"))
match <- as.tibble(match)
#Add tally to match and then collapse by unique vdj
match<- match %>% filter(VJ_HC.x == VJ_HC.y) #keep matching VJ
conv<- match %>% arrange(desc(clonal_abundance)) %>% select(ID ,Source,clonal_abundance, c_gene, name, VJ_HC.x, VJ_HC.y,total_mut, junction_aa.x, junction_aa.y, epitope, neutralising, Sources)  %>% as.data.frame()

write.csv(conv, "/well/combat/users/vkh192/repertoire/out/convergent_sars.csv")

dist <- stringdist::stringdistmatrix(conv$junction_aa.x, conv$junction_aa.y, method="lv")
rownames(dist) <- make.unique(as.character(conv$Source))
colnames(dist) <- make.unique(as.character(conv$junction_aa.y))

rowanno <- conv %>% select(neutralising) %>% mutate(neutralising= ifelse(is.na(.), "non-neut", neutralising))
rownames(rowanno) <- rownames(dist)

colanno <- conv %>% select(c_gene, total_mut, name, Sources) 
rownames(colanno) <- colnames(dist)

png("/well/combat/users/vkh192/repertoire/out/cov_heatmap.png",type=c("cairo"), height = 10, width = 20, units = 'in', res = 400)
pheatmap::pheatmap(dist, annotation_row = rowanno, annotation_col=colanno)
dev.off()

nrow(rowanno)

nrow(dist)
length(conv$neutralising)

#by AA
  
combat <- Heavy %>% select(barcode_seq, sequence_aa, Source, COMBAT_ID) %>%
  rename(ID=barcode_seq, Hchain_aa=sequence_aa, antigen_name=Source, antigen_species=COMBAT_ID)
combat$antigen_species <- "COMBAT"
combine <- rbind(combat,sarscov2)
combine$num <- c(1:nrow(combine))
combine$Hchain_aa <- str_replace_all(combine$Hchain_aa, "[[:punct:]]", "")
#cdhit function, calling on cdhit functions in src
cdhit <- function (seqs, options, name = "CD-Hit", showProgress = interactive()) 
{
  seqs <- Biostrings::AAStringSet(seqs)
  options$i <- tempfile()
  writeXStringSet(seqs, options$i)
  on.exit(unlink(options$i))
  switch(class(seqs), AAStringSet = BcellR:::cdhitC(options, name, showProgress) + 
           1, DNAStringSet = BcellR:::cdhitestC(options, name, showProgress) + 
           1, stop("seqs must be either AAStringSet or DNAStringSet"))
}

#concatenate V and J segments into a new vector

vjseq <- Biostrings::AAStringSet(combine$Hchain_aa)
length(vjseq)
names(vjseq) <- combine$num
vjseq@seq %>% summary
#Set CDhit opts, for now not changeable to user directly. Can incorporate if needed
cdhitOpts <- setNames(as.list(c(0.97, 0, 1)), c("c", "M", "AL"))
cdhitOpts <- lapply(cdhitOpts, as.character)

#run cdhit to cluster vj segments of 80% identity
cdhitclusters<- cdhit(vjseq, cdhitOpts, 'Preclustering')

#extract corresponding whole vdj
vdj <- data.table::as.data.table(cbind(names(vjseq), unlist(cdhitclusters), combine$ID, combine$Hchain_aa, as.character(combine$antigen_name), combine$antigen_species))
colnames(vdj) <- c("sequence_number", "cdhit_cluster", "ID", "Hchain_aa", "antigen_name", "antigen_species")

vdj %>% group_by(antigen_species, cdhit_cluster) %>% add_tally %>% arrange(desc(n)) %>% as.data.frame %>% head

vdj %>% 
  group_by(cdhit_cluster) %>% 
  filter(dplyr::n_distinct(antigen_species) == 2) %>%
  arrange(cdhit_cluster) %>% as.data.frame %>% head(200)




