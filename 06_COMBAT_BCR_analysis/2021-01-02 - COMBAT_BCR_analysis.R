##Analysis/plotting of COMBAT BCR single and bulk repertoires
##Author: Bo Sun
##Group: Bashford-Rogers
##Department: Wellcome Center of Human Genetics
##Institute: Universtit of Oxford
##Contact: bo.sun@well.ox.ac.uk

##libraries----
library(tidyverse)
library(ggplot2); theme_set(theme_classic()+
                              theme(axis.ticks.length=unit(.15, "cm"), 
                                    axis.text.x = element_text(angle=45, hjust=1), 
                                    axis.line = element_line(colour = 'black', size = 0.2), 
                                    axis.ticks = element_line(color="black", size=0.2), strip.background =element_rect(fill="black"),
                                      strip.text = element_text(colour = 'white')))
library(RColorBrewer)
library(mosaic)
library(grid)
library(ggthemes)
library(shadowtext)
library(ggrepel)
library(rstatix)
library(alakazam)
library(purrr)
library(MKmisc)
library(igraph)
library(ggnetwork)
library(ITNr)
library(intergraph)
library(ggnet2)
library(ggnet)
library(GGally)
library(UpSetR)
library(ggnewscale)
library(stringdist)
library(ggmsa)
library(msa)

#load Bcell_df
Bcell_df <- data.table::fread("/Users/bosun/Documents/GitHub/COMBAT_rep/data/Bcell_df_v001.csv")
Light <- data.table::fread("/Users/bosun/Documents/GitHub/COMBAT_rep/data/Bcell_df_LC_v001.csv")

##Set plotting params ----

group.colors <- c(CC = "#D55E00", CS = "#E69F00", CM ="#F0E442", CComm = "#56B4E9", CConv = "#0072B2", HV = "#009E73", Flu = "#999999", Sepsis = "#CC79A7", LCC = "#000000")

nb.cols <- 27
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)
theme_Publication <- function(base_size=14, base_family="helvetica") {
  
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}
# Function to produce summary statistics (mean and +/- sd)----
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#get samples that pass QC----
samplenames <- Bcell_df %>% 
  group_by(COMBAT_ID_Time, Source) %>%
  add_tally(name="n_bcr") %>%
  filter(n_bcr>30)%>%
  ungroup() %>%
  filter(!grepl("LDN", Source)) %>%
  group_by(COMBAT_ID) %>% 
  dplyr::slice(which.min(Timepoint)) %>%
  ungroup() %>%
  select(COMBAT_ID_Time)

##Further clean of data and set factors----
Heavy <- Bcell_df[c(Bcell_df$COMBAT_ID_Time %in% samplenames$COMBAT_ID_Time), ]
Heavy <- Heavy %>% 
  filter(!grepl("Flu",Source_abrev)) #Remove flu cohort from analysis due to poor B cell recovery
  
source <- Heavy %>% 
  ungroup() %>% 
  group_by(COMBAT_ID_Time, RNASeq_sample_ID,Source_abrev) %>% 
  add_tally() %>% 
  ungroup() %>% 
  select(COMBAT_ID_Time, RNASeq_sample_ID,Source_abrev,Age, n, sampled_at_max_severity) %>% 
  unique() %>% 
  as.data.frame() %>%
  arrange(desc(n))
source2 <- Heavy %>% 
  select(Source_abrev,Source) %>% 
  unique() %>% 
  as.data.frame()

source$Source_abrev <- factor(source$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))
source2$Source_abrev <- factor(source2$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))
Heavy$Source_abrev <- factor(Heavy$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))
Heavy$pseudobulk <- factor(Heavy$pseudobulk, levels = c("B.NAIVE", "B.INT","B.MEM", "PB"))

##Sanity checks----
Heavy %>% 
  filter(!is.na(pseudobulk)) %>%
  group_by(COMBAT_ID_Time, pseudobulk) %>% 
  add_tally() %>%
  select(COMBAT_ID_Time, pseudobulk, n) %>%
  unique %>%
  as.data.frame()%>%
  arrange(n)

Heavy %>% 
  group_by(COMBAT_ID_Time, pseudobulk) %>% 
  add_tally() %>% 
  select(n) %>% 
  select(n) %>%
  unique %>% 
  arrange(n) %>%
  as.data.frame()%>%
  head(50)

##Generate umap figure----
##General umap plot

#Merge switched and unswitched populations for visualisation
Heavy <- Heavy %>%
  mutate(new_anno = ifelse(grepl("B.NAIVE|B.NAIVE.IgDlo|B.NAIVE.CD1c|B.TRANSIT.CD10|B.NAIVE.IFN.resp", name.unique), "Näive B",
                            ifelse(grepl("B.SW.MEM.IFN.resp|B.SW.MEM", name.unique), "Switched memory B", 
                                  ifelse(grepl("B.UNSW.MEM|B.int.2.unsw", name.unique), "Unswitched memory B", 
                                         ifelse(grepl("PB.IFN.resp|PB.mitohi|PB", name.unique), "Plasmablast", 
                                                ifelse(grepl("PB.cyc|B.cyc", name.unique), "Plasmablast", 
                                                        ifelse(grepl("B.int.2.early.act/sw|B.int.1.early.act/sw|B.mitohi|B.int.1.IFN.resp|B.int.1.early.act|B.int.2.IFN.resp|B.int.2.early.act.IFN.resp", name.unique), "Intermediate B", NA)))))))

Heavy$new_anno <- factor(Heavy$new_anno, levels=c("Näive B", "Intermediate B", "Switched memory B", "Unswitched memory B", "Plasmablast"))

umap_plot<- Heavy %>% 
  ggplot(aes(UMAP_1, UMAP_2))+
  geom_point(aes(color=new_anno), size=1, shape=16, alpha=0.5, lwd=0.5)+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggsci::scale_color_jama()+
  ggsci::scale_fill_jama()+
  theme_void()+
  theme(legend.position = "none")

pdf("COMBAT_figure_plots/umap_plot.pdf", height = 8, width = 10)
print(umap_plot)
dev.off()

#umap of umis 
umap_umi_plot<- Heavy %>% 
  ggplot(aes(UMAP_1, UMAP_2))+
  geom_point(aes(color= log(umis)), size=1, shape=16, alpha=0.8, lwd=0.5)+
  viridis::scale_color_viridis()+
  theme_void()

pdf("COMBAT_figure_plots/umap_umi_plot.pdf", height = 8, width = 10)
print(umap_umi_plot)
dev.off()

#plot proportions    
clusterprop <- Heavy %>%
  group_by(COMBAT_ID_Time)%>%
  add_tally(name = "sum") %>%
  group_by(COMBAT_ID_Time, pseudobulk) %>% 
  add_tally(name="sum_p") %>%
  select(COMBAT_ID_Time, Source_abrev, pseudobulk,sum,sum_p) %>% 
  unique %>%
  mutate(prop=sum_p/sum) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  select(-sum, -sum_p) %>%
  spread(pseudobulk, value=prop)

clusterprop[is.na(clusterprop)] <- 0

clusterprop <- clusterprop %>%
  reshape2::melt(id.vars=c("COMBAT_ID_Time", "Source_abrev"))


clusterprop <- clusterprop %>%
  mutate(variable= ifelse(grepl("NAIVE", variable), "Näive B", ifelse(grepl("PB", variable), "Plasmablast", ifelse(grepl("INT", variable), "Intermediate B", "Memory B"))))
clusterprop$variable <- factor(clusterprop$variable, levels=c("Näive B", "Intermediate B","Memory B",  "Plasmablast"))
cluspropplot <- clusterprop %>%
  ggplot(aes(Source_abrev, value))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.6)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  facet_wrap(~variable, nrow=1)+
  xlab("")+
  ylab("Propotion of all B cells")

clustpropdunn<- clusterprop %>%
  group_by(variable) %>%
  rstatix::dunn_test(value~Source_abrev, p.adjust.method = "fdr")
clustpropdunn$variable <- factor(clustpropdunn$variable, levels=c("Näive B", "Intermediate B","Memory B",  "Plasmablast"))

cluspropplot <- cluspropplot+ ggpubr::stat_pvalue_manual(clustpropdunn %>% filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 0.7, step.increase = 0.1,)

pdf("COMBAT_figure_plots/cluster_proportions.pdf", height = 4, width = 9)
print(cluspropplot)
dev.off()

readr::write_csv(clusterprop, "out/cluster_prop.csv")
readr::write_csv(clustpropdunn, "out/cluster_prop_dunn.csv")

##Alpha diversity scores----
sample_curve <- alakazam::alphaDiversity(Heavy, group="Source_abrev", clone="clone_per_replicate",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=200)

adivplot <- alakazam::plot(sample_curve, colors=group.colors, main_title="Alpha diversity", 
     legend_title="")+
  theme_classic()+
  theme(axis.ticks.length=unit(.15, "cm"), 
        axis.text.x = element_text(angle=45, hjust=1), 
        axis.line = element_line(colour = 'black', size = 0.2), 
        axis.ticks = element_line(color="black", size=0.2), strip.background =element_rect(fill="black"),
        strip.text = element_text(colour = 'white'))

pdf("COMBAT_figure_plots/alpha_diversity.pdf", height = 4, width = 4)
print(adivplot)
dev.off()

saveRDS(sample_curve, "out/alpha_diversity_samplecurve.rds")

##Custom diversity metrics----
Heavy %>% 
  group_by(COMBAT_ID_Time, pseudobulk) %>% 
  add_tally() %>%
  select(COMBAT_ID_Time,Source, pseudobulk, n) %>%
  unique %>%
  as.data.frame() %>%
  arrange(n)

sampletest <- Heavy %>% 
  group_by(COMBAT_ID_Time, pseudobulk) %>% 
  add_tally() %>% 
  filter(n>=10) 

samplepass <- sampletest %>%
  select(COMBAT_ID_Time, Source, pseudobulk) %>%
  unique() %>%
  group_by(COMBAT_ID_Time) %>%
  add_tally(name="pseudo_bulk_count") %>%
  select(COMBAT_ID_Time, pseudobulk, Source, pseudo_bulk_count) %>%
  unique() %>%
  filter(pseudo_bulk_count == 4) %>%
  select(COMBAT_ID_Time) %>%
  unique
  
diversity_samples <- Heavy[c(Heavy$COMBAT_ID_Time %in% samplepass$COMBAT_ID_Time),]
namelist <- diversity_samples %>% 
  group_by(COMBAT_ID_Time, pseudobulk) %>% 
  add_tally() %>% 
  filter(n>=10) %>% 
  group_split() 


source("Functions/Get_diversity.R")
#sanity check filter 
lengths <- lapply(namelist, function(x){x<-nrow(x)})
lengths<-unlist(lengths)
min(lengths)

diversity <- get_bcr_metrics(namelist, method="all", repertoire_id = "COMBAT_ID_Time", cell_subset = "pseudobulk", clone_id = "clone_per_replicate", subsample=8)

saveRDS(diversity, "out/diversity.rds")

##Diversity stats----
boot.CDI.kw <- diversity[[3]] %>% 
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  group_by(pseudobulk) %>%
  rstatix::kruskal_test(mean.CDI ~ Source_abrev)

boot.CDI.ancova.age <- diversity[[3]] %>% 
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  group_by(pseudobulk) %>%
  rstatix::anova_test(mean.CDI ~ Source_abrev, covariate=Age) 

boot.shannon.dunn <- diversity[[1]] %>% 
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  group_by(pseudobulk) %>%
  rstatix::dunn_test(mean.shannon ~ Source_abrev, p.adjust.method = "fdr") %>%
  rstatix::add_x_position(x = "Source_abrev")

boot.CEI.dunn <- diversity[[2]] %>% 
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  group_by(pseudobulk) %>%
  rstatix::dunn_test(mean.CEI ~ Source_abrev, p.adjust.method = "fdr") %>%
  rstatix::add_x_position(x = "Source_abrev")

boot.CDI.dunn <- diversity[[3]] %>% 
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  group_by(pseudobulk) %>%
  rstatix::dunn_test(mean.CDI ~ Source_abrev, p.adjust.method = "fdr") %>%
  rstatix::add_x_position(x = "Source_abrev")

boot.shannon.dunn %>%
  filter(p.adj.signif != "ns")
boot.CEI.dunn %>%
  filter(p.adj.signif != "ns")
boot.CDI.dunn %>%
  filter(p.adj.signif != "ns")

readr::write_csv(boot.CDI.kw, "out/CDI.kw.csv")
readr::write_csv(boot.CDI.ancova.age, "out/CDI.ancova.age.csv")
readr::write_csv(boot.CDI.dunn, "out/CDI.dunn.csv")

##Here we can see significance is only reached in Plasmablast population.. will plot ----
cdi_boxplot <- diversity[[3]] %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  filter(pseudobulk == "PB") %>%
  mutate(labels="Plasmablast") %>%
  mutate(group_1=Source_abrev) %>%
  ggplot(aes(Source_abrev, mean.CDI))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.6)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  xlab("")+
  ylab("Mean CDI")+
  facet_wrap(~labels)+
  theme(legend.position = "none")

cdi_boxplot<- cdi_boxplot + ggpubr::stat_pvalue_manual(boot.CDI.dunn %>%   filter(pseudobulk == "PB") %>%
filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 0.35, step.increase = 0.06,)
cdi_boxplot

pdf("COMBAT_figure_plots/cdi_plasmablasts.pdf", height = 4, width=2.5)
print(cdi_boxplot)
dev.off()

shannon_boxplot <- diversity[[1]] %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  filter(pseudobulk == "PB") %>%
  mutate(labels="Plasmablast") %>%
  mutate(group_1=Source_abrev) %>%
  ggplot(aes(Source_abrev, mean.shannon))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.6)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  xlab("")+
  ylab("Mean Shannon")+
  facet_wrap(~labels)+
  theme(legend.position = "none")

shannon_boxplot<- shannon_boxplot + ggpubr::stat_pvalue_manual(boot.shannon.dunn %>%   filter(pseudobulk == "PB") %>%
                                                         filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 2.2, step.increase = 0.06,)
shannon_boxplot

pdf("COMBAT_figure_plots/shannon_plasmablasts.pdf", height = 4, width=2.5)
print(shannon_boxplot)
dev.off()

cei_boxplot <- diversity[[2]] %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  filter(pseudobulk == "PB") %>%
  mutate(labels="Plasmablast") %>%
  mutate(group_1=Source_abrev) %>%
  ggplot(aes(Source_abrev, mean.CEI))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.2, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.6)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  xlab("")+
  ylab("Mean CEI")+
  facet_wrap(~labels)+
  theme(legend.position = "none")

cei_boxplot<- cei_boxplot + ggpubr::stat_pvalue_manual(boot.CEI.dunn %>%   filter(pseudobulk == "PB") %>%
                                                         filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 0.25, step.increase = 0.06,)
cei_boxplot

pdf("COMBAT_figure_plots/CEI_plasmablasts.pdf", height = 4, width=2.5)
print(CEI_boxplot)
dev.off()

##plotname cdi_cei_plot x2 horizontal
pdf("COMBAT_figure_plots/cei_cdi_plasmablasts.pdf", height = 4, width=4.5)
patchwork::wrap_plots(CEI_boxplot|cdi_boxplot, nrow=1)
dev.off()

##Set up individual statistics----
#c_gene----
c_gene_prop <- Heavy %>%
  filter(c_gene != "None") %>%
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum_c_gene") %>%
  group_by(COMBAT_ID_Time, pseudobulk, c_gene) %>% 
  add_tally(name="sum_isotype") %>%
  select(COMBAT_ID_Time, pseudobulk, c_gene,sum_c_gene,sum_isotype) %>% 
  unique %>%
  mutate(prop_isotype=sum_isotype/sum_c_gene) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  select(-sum_c_gene, -sum_isotype) %>%
  spread(c_gene, value= prop_isotype) 

c_gene_prop[is.na(c_gene_prop)] <- 0

c_gene_prop <- c_gene_prop %>%
  reshape2::melt(., id=c("COMBAT_ID_Time", "pseudobulk"))

readr::write_csv(c_gene_prop, "out/c_gene_prop.csv")

#c_gene stats----
c_genelist <- c_gene_prop%>% 
  dplyr::left_join(source, by=c("COMBAT_ID_Time"="COMBAT_ID_Time")) %>%
  filter(!(grepl("Flu", Source_abrev))) %>%
  group_by(pseudobulk, variable) %>%
  group_split 

#c_gene kw screen
cgene.kw <- lapply(c_genelist, function(y){
  if(sum(y$value) != 0){
    kw<- rstatix::kruskal_test(y, value ~ Source_abrev) 
    return(kw)  
  }else{
    NULL
  }
})

#non significant 

#No significance found per cell subset----
cgene_plot<- c_gene_prop %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  filter(variable != "No_c_gene") %>%
  filter(!grepl("IGHE", variable)) %>%
  filter(!grepl("Flu",Source_abrev)) %>%
  ggplot(aes(Source_abrev, value, fill=Source_abrev)) +
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.1, shape=21)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
  geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  facet_grid(vars(pseudobulk), vars(variable), scales = "free_y")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  xlab("")+
  ylab("Proportion of repertoire")
cgene_plot

##c_gene plot 
pdf("COMBAT_figure_plots/c_gene_supplementary.pdf", height=10, width=10)
cgene_plot
dev.off()

#Proportion clonal ----
clonal_prop <- Heavy %>% 
  mutate(clonal=ifelse(clonal_abundance>1, "Clonal", "Non-clonal")) %>%
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  group_by(COMBAT_ID_Time, pseudobulk, clonal) %>% 
  add_tally(name="sum_clonal") %>%
  select(COMBAT_ID_Time, pseudobulk,Source_abrev, clonal,sum,sum_clonal) %>% 
  unique %>%
  mutate(prop_clonal=sum_clonal/sum) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  select(-sum, -sum_clonal) %>%
  filter(!grepl("Non", clonal)) %>%
  ungroup %>%
  unique  %>%
  spread(pseudobulk, value= prop_clonal) 

clonal_prop[is.na(clonal_prop)] <- 0

clonal_prop <- clonal_prop %>%
  reshape2::melt(., id=c("COMBAT_ID_Time", "Source_abrev","clonal"))

readr::write_csv(clonal_prop, "out/clonal_morethan2_prop.csv")

##Proportion clone >2 stats ----
clonalkw <- clonal_prop %>%
  group_by(variable) %>%
  rstatix::kruskal_test(value ~Source_abrev) 
#Significance in PB
clonaldunn <- clonal_prop %>%
  group_by(variable) %>%
  rstatix::dunn_test(value ~Source_abrev, p.adjust.method = "fdr") %>%
  rstatix::add_x_position(x="Source_abrev")

#Clonal >2 plot
clonal_prop_boxplot <- clonal_prop %>%
  filter(grepl("PB", variable)) %>%
  ggplot(aes(Source_abrev, value))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.6)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  facet_wrap(~variable)+
  xlab("")+
  ylab("Proportion of repertoire (clonesize >=2)")+
  theme(legend.position = "none")
clonal_prop_boxplot <- clonal_prop_boxplot + ggpubr::stat_pvalue_manual(clonaldunn %>% filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 0.82, step.increase = 0.1,)

pdf("COMBAT_figure_plots/clonal_prop_pb.pdf", height=4, width=2.5)
clonal_prop_boxplot
dev.off()

pdf("COMBAT_figure_plots/clonal_prop_pb_cdi_cei.pdf", height=4, width=7)
patchwork::wrap_plots(clonal_prop_boxplot|CEI_boxplot|cdi_boxplot, nrow=1)
dev.off()

##prop clonal c_gene
clonal_prop_c <- Heavy %>% 
  filter(clonal_abundance>1)%>%
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  group_by(COMBAT_ID_Time, pseudobulk, c_gene) %>% 
  add_tally(name="sum_clonal") %>%
  select(COMBAT_ID_Time,Source_abrev,sum,sum_clonal, c_gene, clonal_abundance, pseudobulk) %>% 
  unique %>%
  mutate(prop_clonal=sum_clonal/sum) %>%
  select(-sum, -sum_clonal) %>%
  spread(c_gene, value= prop_clonal) 

clonal_prop_c[is.na(clonal_prop_c)] <- 0

clonal_prop_c <- clonal_prop_c %>%
  reshape2::melt(., id=c("COMBAT_ID_Time", "Source_abrev", "pseudobulk", "clonal_abundance"))

#clonal_c_gene stats----
clonal_c_kw<- clonal_prop_c %>% 
  group_by(pseudobulk, variable) %>%
  rstatix::kruskal_test(value~Source_abrev) %>%
  filter(p < 0.05)

#c_gene dunn screen
clonal_cgene.dunn <- lapply(c_genelist, function(y){
  if(sum(y$value) != 0){
    dunn<- rstatix::dunn_test(y, value ~ Source_abrev, p.adjust.method = "fdr") 
    dunn <- dunn %>% mutate(variable= paste(unique(y$pseudobulk), unique(y$variable), sep="_"))
    return(dunn)  
  }else{
    NULL
  }
})

clonal_prop_c %>%
  filter(grepl("PB", pseudobulk)) %>%
  ggplot(aes(Source_abrev, value))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~variable)

##Dunn does not show pairwise significance

#IGHV----
#create proportion table of Heavy V genes. Here each V gene usage is represented as the proportion of BCRs using that V gene in an individuals whole repertoire
ighv <- Heavy %>%
  filter(!grepl("Flu", Source)) %>%
  mutate(Source_global = ifelse(grepl("COVID", Source), "COVID", ifelse(grepl("HV", Source), "HV", "Sepsis")))

HeavyV <- as.data.frame.matrix(table(ighv$COMBAT_ID_Time, ighv$v_gene_HC) %>% 
                                 prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.)) %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>% 
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID)%>% 
  reshape2::melt(., id=c("COMBAT_ID_Time","Source_abrev", "Age")) 
#Reorder table for kruskal wallis test
heavy_v_list <- HeavyV %>% 
  group_by(variable) %>%
  group_split() %>%
setNames(unique(HeavyV$variable))

#CDI kw+dunn
igh_kw <- lapply(heavy_v_list, function(x){
  x %>% 
    kruskal_test(value~Source_abrev)
})

#Post hoc Dunn
igh_dunn <- lapply(heavy_v_list, function(x){
  x %>%
    rstatix::dunn_test(value~Source_abrev, p.adjust.method = "fdr", detailed = FALSE) %>%
    rstatix::add_xy_position(x = "Source_abrev") %>%
    mutate(.y. = unique(x$variable))
})

heavy_v_means_list <- HeavyV %>%
  group_by(variable) %>%
  group_split() %>%
setNames(unique(HeavyV$variable))

long_mean_logfc <- lapply(heavy_v_means_list, function(x){
  pairwise.fc(x$value, x$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
  data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
  mutate(.y.=unique(x$variable)) %>%
    mutate(condition=paste(condition1, .y., sep="-")) %>%
  separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
})

long_mean_logfc <- data.table::rbindlist(long_mean_logfc)
igh_dunn_bound <- data.table::rbindlist(igh_dunn)

igh_dunn_bound <- igh_dunn_bound %>% 
  left_join(long_mean_logfc, by=c("group1", "group2", ".y.")) 

ighp1 <- igh_dunn_bound %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.6)+
  geom_point(shape=16, alpha=0.8)+
  ggrepel::geom_text_repel(data=filter(igh_dunn_bound, p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  theme(legend.position = "none")+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))

readr::write_csv(igh_dunn_bound %>% select(-groups), "out/igh_dunn_single_cell.csv")

#IGHL----

#create proportion table of Heavy V genes. Here each V gene usage is represented as the proportion of BCRs using that V gene in an individuals whole repertoire
ighl_group <- Light %>%
  filter(!grepl("Flu", Source)) %>%
  mutate(Source_global = ifelse(grepl("COVID", Source), "COVID", ifelse(grepl("HV", Source), "HV", "Sepsis"))) %>%
  select(COMBAT_ID_Time, Source_global)
ighl <- Light %>%
  filter(!grepl("Flu", Source)) %>%
  mutate(Source_global = ifelse(grepl("COVID", Source), "COVID", ifelse(grepl("HV", Source), "HV", "Sepsis")))

LightV <- as.data.frame.matrix(table(ighl$COMBAT_ID_Time, ighl$v_gene_LC) %>% 
                                 prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.)) %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>% 
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID)%>% 
  reshape2::melt(., id=c("COMBAT_ID_Time","Source_abrev", "Age")) 
#Reorder table for kruskal wallis test

light_v_list <- LightV %>% 
  group_by(variable) %>%
  group_split() %>%
  setNames(unique(LightV$variable))

#CDI kw+dunn

igl_kw <- lapply(light_v_list, function(x){
  x %>% 
    kruskal_test(value~Source_abrev)
})
igl_ancov <- lapply(light_v_list, function(x){
  x %>% 
    anova_test(value~Source_abrev, covariate=Age)
})

#Post hoc Dunn

igl_dunn <- lapply(light_v_list, function(x){
  x %>%
    rstatix::dunn_test(value~Source_abrev, p.adjust.method = "fdr", detailed = FALSE) %>%
    rstatix::add_xy_position(x = "Source_abrev") %>%
    mutate(.y. = unique(x$variable))
})

light_v_means_list <- LightV %>%
  group_by(variable) %>%
  group_split() %>%
  setNames(unique(LightV$variable))

long_igl_mean_logfc <- lapply(light_v_means_list, function(x){
  pairwise.fc(x$value, x$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
    data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
    mutate(.y.=unique(x$variable)) %>%
    mutate(condition=paste(condition1, .y., sep="-")) %>%
    separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
})

long_igl_mean_logfc <- data.table::rbindlist(long_igl_mean_logfc)
igl_dunn_bound <- data.table::rbindlist(igl_dunn)

igl_dunn_bound <- igl_dunn_bound %>% 
  left_join(long_igl_mean_logfc, by=c("group1", "group2", ".y.")) 
iglp1 <- igl_dunn_bound %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.6)+
  geom_point(shape=16, alpha=0.8)+
  ggrepel::geom_text_repel(data=filter(igl_dunn_bound, p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  guides(color=guide_legend(title="FDR"))+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))

readr::write_csv(igl_dunn_bound %>% select(-groups), "out/igl_dunn_single_cell.csv")

pdf("COMBAT_figure_plots/combined_ighv_igl_volcano.pdf", height=6, width=11)
patchwork::wrap_plots(ighp1|iglp1, nrow=1)
dev.off()


##alternative plots
igh_dunn_bound %>%
  mutate(condition=gsub("-.*", "", condition)) %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(x=.y., y=log2(fold_change), color=col))+
  geom_point(shape=19, color="black", alpha=0.6)+
  geom_point(shape=16, alpha=0.8)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  xlab("avg_logFC of repertoire proportion")+
  guides(color=guide_legend(title="FDR"))+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))+
  coord_flip()+
  facet_wrap(~condition)

#IGHVJ----
#create proportion table of Heavy V genes. Here each V gene usage is represented as the proportion of BCRs using that V gene in an individuals whole repertoire
ighvj <- Heavy %>%
  mutate(vj_pair= paste(v_gene_HC, j_gene_HC, sep="_")) 

ighvj <- ighvj %>%
  filter(!grepl("Flu", Source)) %>%
  mutate(Source_global = ifelse(grepl("COVID", Source), "COVID", ifelse(grepl("HV", Source), "HV", "Sepsis")))

ighVJ <- as.data.frame.matrix(table(ighvj$COMBAT_ID_Time, ighvj$vj_pair) %>% 
                                 prop.table(margin = 1)) %>% 
  mutate(COMBAT_ID_Time= rownames(.)) %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>% 
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID)%>% 
  reshape2::melt(., id=c("COMBAT_ID_Time","Source_abrev", "Age")) 

#Reorder table for kruskal wallis test
vj_v_list <- ighVJ %>% 
  group_by(variable) %>%
  group_split() %>%
  setNames(unique(ighVJ$variable))

#CDI kw+dunn
vj_kw <- lapply(vj_v_list, function(x){
  x %>% 
    kruskal_test(value~Source_abrev)
})
vj_ancov <- lapply(vj_v_list, function(x){
  x %>% 
    anova_test(value~Source_abrev, covariate=Age)
})

#Post hoc Dunn
igvj_dunn <- lapply(vj_v_list, function(x){
  x %>%
    rstatix::dunn_test(value~Source_abrev, p.adjust.method = "fdr", detailed = FALSE) %>%
    rstatix::add_xy_position(x = "Source_abrev") %>%
    mutate(.y. = unique(x$variable))
})

igvj_dunn %>%
  data.table::rbindlist() %>%
  filter(!grepl("ns", p.adj.signif))
igvj_dunn2 %>%
  data.table::rbindlist() %>%
  filter(!grepl("ns", p.adj.signif))

igvj_dunn<- ighVJ %>% 
  group_by(variable) %>%
  rstatix::dunn_test(value~Source_abrev, p.adjust.method = "fdr", detailed = FALSE)

hvj_means_list <- ighVJ %>%
  group_by(variable) %>%
  group_split() %>%
  setNames(unique(ighVJ$variable))

long_hvj_mean_logfc <- lapply(hvj_means_list, function(x){
  pairwise.fc(x$value, x$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
    data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
    mutate(.y.=unique(x$variable)) %>%
    mutate(condition=paste(condition1, .y., sep="-")) %>%
    separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
})

long_hvj_mean_logfc <- data.table::rbindlist(long_hvj_mean_logfc)
igvj_dunn_bound <- data.table::rbindlist(igvj_dunn)

igvj_dunn_bound <- igvj_dunn %>% 
  left_join(long_hvj_mean_logfc, by=c("group1", "group2", "variable"=".y.")) 
igvjp1 <- igvj_dunn_bound %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.6)+
  geom_point(shape=16, alpha=0.8)+
  ggrepel::geom_text_repel(data=filter(igvj_dunn_bound, p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3, 
                           # Repel away from the left edge, not from the right.
                           xlim = c(NA, Inf),
                           # Do not repel from top or bottom edges.
                           ylim = c(-Inf, Inf))+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  guides(color=guide_legend(title="FDR"))+
  theme(legend.position = "none")

vj_dotplot<- igvj_dunn_bound %>%
  mutate(group1vs=paste(group1, "vs.", sep=" ")) %>%
  filter(!grepl("ns", p.adj.signif)) %>%
  filter(!grepl("Inf", fold_change)) %>%
  filter(fold_change != 0) %>%
  ggplot(aes(variable, group2))+
  geom_point(aes(size=log2(fold_change),color=log2(fold_change)), alpha=0.8, shape=19)+
  geom_text(aes(label=p.adj.signif), vjust=0.8)+
  facet_wrap(~group1vs, nrow=1)+
  ggsci::scale_color_gsea()+
  coord_flip()+
  ylab("")+
  xlab("IGHV-J pairing")

readr::write_csv(igvj_dunn_bound, "out/ighvj_dunn_single_cell.csv")

pdf("COMBAT_figure_plots/vj_dotplot.pdf", height=4, width=8)
vj_dotplot
dev.off()
#Median SMH per pseudobulk ----
tmut <- data.table::setDT(Heavy)[,list(Mean.total.mut=mean(total_mut), Max.total.mut=max(total_mut), Min.total.mut=min(total_mut), Median.total.mut=as.numeric(median(total_mut)), Std.total.mut=sd(total_mut)), by=c("COMBAT_ID_Time", "pseudobulk")]

tmut %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time")) %>%
  ggplot(aes(Source_abrev, Median.total.mut, fill=Source_abrev))+
  geom_boxplot(alpha=0.5, size=0.2, outlier.shape = NA)+
  geom_jitter(width=0.1,shape=21, size=1)+
  facet_wrap(~pseudobulk)+
  scale_fill_manual(values=group.colors)

#Reorder table for kruskal wallis test
tmut2 <- tmut %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>% 
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID)
tmut_list <- tmut2 %>%
  group_by(pseudobulk) %>%
  group_split() 

#CDI kw+dunn
tmut_kw <- lapply(tmut_list, function(x){
  x %>% 
    kruskal_test(Median.total.mut~Source_abrev)
})

#Post hoc Dunn
tmut_dunn <-tmut %>%
  group_by(pseudobulk) %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  filter(Median.total.mut != "B.NAIVE") %>%
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID) %>%
  rstatix::dunn_test(Median.total.mut~Source_abrev, p.adjust.method = "fdr", detailed = FALSE)%>%
  rstatix::add_xy_position(x = "Source_abrev")

readr::write_csv(tmut_dunn %>% select(-groups), "out/Median_mut_dunn_sc.csv")

tmut_plot<- tmut %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  filter(pseudobulk != "B.NAIVE") %>%
  ggplot(aes(Source_abrev, Median.total.mut))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, lwd=0) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.5)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  facet_wrap(~pseudobulk, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  xlab("")+
  ylab("Median total IGHV mutations")
tmut_plot + ggpubr::stat_pvalue_manual(tmut_dunn %>% filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 30, step.increase = 0.07)

tmut_dunn %>% filter(grepl("PB", pseudobulk))

#SMH <4 prop ---
lowmut_clone <- Heavy %>% 
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  mutate(lowmut= total_mut <4)%>%
  group_by(COMBAT_ID_Time, pseudobulk, lowmut) %>% 
  add_tally(name="low_mut") %>%
  ungroup %>%
  select(COMBAT_ID_Time,RNASeq_sample_ID, pseudobulk, lowmut ,sum,low_mut, Source_abrev) %>% 
  unique %>%
  mutate(prop_ighv=low_mut/sum) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  select(-sum, -low_mut) %>%
  spread(lowmut, value=prop_ighv) %>%
  dplyr::rename(high_mut=`FALSE`, low_mut=`TRUE`) 

lowmut_clone[is.na(lowmut_clone)] <- 0
lowmut_clone %>%
  ggplot(aes(Source_abrev, low_mut, fill=Source_abrev))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(shape=21, size=0.9, width=0.2)+
  facet_wrap(~pseudobulk, scales = "free_y")+
  scale_fill_manual(values=group.colors)

#stats
#Post hoc Dunn
lowmut_dunn <-lowmut_clone %>%
  group_by(pseudobulk) %>%
  filter(pseudobulk != "B.NAIVE") %>%
  rstatix::dunn_test(low_mut~Source_abrev, p.adjust.method = "fdr", detailed = FALSE)%>%
  rstatix::add_xy_position(x = "Source_abrev")

readr::write_csv(lowmut_dunn, "out/low_mut_clone_prop_dunn.csv")
lowmut_dunn %>% 
  filter(!grepl("ns", p.adj.signif))
#No significance

#Median SMH FREQUENCY per pseudobulk ----
Heavy <- Heavy %>% 
  mutate(smh_freq = total_mut/nchar(v_region))

Heavy$smh_freq
mutfreq <- data.table::setDT(Heavy)[,list(Mean.mut.freq=mean(smh_freq), Max.mut.freq=max(smh_freq), Min.mut.freq=min(smh_freq), Median.mut.freq=as.numeric(median(smh_freq)), Std.mut.freq=sd(smh_freq)), by=c("COMBAT_ID_Time", "pseudobulk")]

mutfreq %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time")) %>%
  ggplot(aes(Source_abrev, Median.mut.freq, fill=Source_abrev))+
  geom_boxplot(alpha=0.5, size=0.2, outlier.shape = NA)+
  geom_jitter(width=0.1,shape=21, size=1)+
  facet_wrap(~pseudobulk)+
  scale_fill_manual(values=group.colors)
Heavy %>%
  ggplot(aes(smh_freq)) +
  geom_histogram()+
  facet_wrap(~pseudobulk)

#Reorder table for kruskal wallis test
mutfreq %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>% 
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID)%>%
  group_by(pseudobulk) %>%
  kruskal_test(Median.mut.freq~Source_abrev)

mutfreq_dunn <-mutfreq %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID, -Age) %>%
  ungroup %>% 
  filter(Median.mut.freq != 0) %>%
  group_by(pseudobulk) %>%
  group_split()

lapply(mutfreq_dunn, function(x){
  x %>%
    rstatix::dunn_test(Median.mut.freq ~ Source_abrev, p.adjust.method = "fdr", detailed = FALSE) %>%
    rstatix::add_xy_position(x = "Source_abrev") %>%
    filter(p.adj.signif != "ns")
})
  rstatix::dunn_test(Median.mut.freq ~ Source_abrev, p.adjust.method = "fdr", detailed = FALSE) %>%
  rstatix::add_xy_position(x = "Source_abrev")

mutfreq %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  select(-n, -sampled_at_max_severity, -RNASeq_sample_ID, -Age) %>%
  filter(pseudobulk=="PB")%>%
  filter(Median.mut.freq == 0)

mutfreq_plot<- mutfreq %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  filter(pseudobulk != "B.NAIVE") %>%
  filter(Median.mut.freq != 0) %>%
  ggplot(aes(Source_abrev, Median.mut.freq))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, lwd=0) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.5)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  facet_wrap(~pseudobulk, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  xlab("")+
  ylab("Median IGHV mutation frequency")
mutfreq_plot + ggpubr::stat_pvalue_manual(mutfreq_dunn %>% filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 30, step.increase = 0.07)

#SMH <4 prop ---
lowmut_clone <- Heavy %>% 
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  mutate(lowmut= total_mut <4)%>%
  group_by(COMBAT_ID_Time, pseudobulk, lowmut) %>% 
  add_tally(name="low_mut") %>%
  ungroup %>%
  select(COMBAT_ID_Time,RNASeq_sample_ID, pseudobulk, lowmut ,sum,low_mut, Days_symptom_to_sample) %>% 
  unique %>%
  mutate(prop_ighv=low_mut/sum) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  select(-sum, -low_mut) %>%
  spread(lowmut, value=prop_ighv) %>%
  rename(high_mut=`FALSE`, low_mut=`TRUE`)
lowmut_clone %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time")) %>%
  filter(!grepl("Flu", Source_abrev)) %>%
  ggplot(aes(Source_abrev, low_mut, fill=Source_abrev))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(shape=21, size=0.9, width=0.2)+
  facet_wrap(~pseudobulk, scales = "free_y")+
  scale_fill_manual(values=group.colors)

#unmutated class switch
lowmut_clone <- Heavy %>% 
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  mutate(lowmut= total_mut <4)%>%
  group_by(COMBAT_ID_Time, pseudobulk, lowmut) %>% 
  add_tally(name="low_mut") %>%
  ungroup %>%
  select(COMBAT_ID_Time,RNASeq_sample_ID, pseudobulk, lowmut ,sum,low_mut, Days_symptom_to_sample) %>% 
  unique %>%
  mutate(prop_ighv=low_mut/sum) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  select(-sum, -low_mut) %>%
  spread(lowmut, value=prop_ighv) %>%
  rename(high_mut=`FALSE`, low_mut=`TRUE`) %>%
  mutate(low_mut= ifelse(is.na(low_mut), 0, low_mut))
lowmut_clone_cgene <- Heavy %>% 
  group_by(COMBAT_ID_Time, c_gene)%>%
  add_tally(name = "sum") %>%
  mutate(lowmut= total_mut <4)%>%
  group_by(COMBAT_ID_Time, c_gene, lowmut) %>% 
  add_tally(name="low_mut") %>%
  ungroup %>%
  select(COMBAT_ID_Time,RNASeq_sample_ID, c_gene, lowmut ,sum,low_mut, Days_symptom_to_sample) %>% 
  unique %>%
  mutate(prop_ighv=low_mut/sum) %>%
  arrange(c_gene,COMBAT_ID_Time) %>%
  select(-sum, -low_mut) %>%
  spread(lowmut, value=prop_ighv) %>%
  rename(high_mut=`FALSE`, low_mut=`TRUE`) %>%
  mutate(low_mut= ifelse(is.na(low_mut), 0, low_mut))

lowmut_clone_cgene %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time")) %>%
  filter(!grepl("Flu", Source_abrev)) %>%
  ggplot(aes(Source_abrev, low_mut, fill=Source_abrev))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(shape=21, size=0.9, width=0.2)+
  facet_wrap(~c_gene, scales = "free_y")+
  scale_fill_manual(values=group.colors)
##lowmut expanded clone
lowmut_exp_clone <- Heavy %>% 
  filter(c_gene=="IGHG1") %>%
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  mutate(lowmut=total_mut <5)%>%
  group_by(COMBAT_ID_Time, pseudobulk, lowmut) %>% 
  add_tally(name="low_mut") %>%
  ungroup %>%
  select(COMBAT_ID_Time,RNASeq_sample_ID, lowmut ,sum,low_mut, Days_symptom_to_sample, pseudobulk) %>% 
  unique %>%
  mutate(prop_ighv=low_mut/sum) %>%
  arrange(COMBAT_ID_Time) %>%
  select(-low_mut) %>%
  spread(lowmut, value=prop_ighv) %>%
  rename(high_mut=`FALSE`, low_mut=`TRUE`) %>%
  mutate(low_mut= ifelse(is.na(low_mut), 0, low_mut))

lowmut_exp_clone %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time")) %>%
  filter(!grepl("Flu", Source_abrev)) %>%
  ggplot(aes(Source_abrev, low_mut, fill=Source_abrev))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(shape=21, size=0.9, width=0.2)+
  facet_wrap(~pseudobulk, scales = "free_y")+
  scale_fill_manual(values=group.colors)

#stats
#Post hoc Dunn
lowmut_dunn <-lowmut_exp_clone %>%
  group_by(pseudobulk) %>%
  filter(pseudobulk != "B.NAIVE") %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  rstatix::dunn_test(low_mut~Source_abrev, p.adjust.method = "fdr", detailed = FALSE)%>%
  rstatix::add_xy_position(x = "Source_abrev")

lowmut_dunn %>% 
  filter(!grepl("ns", p.adj.signif))


lowmut_exp_clone %>%
  filter(pseudobulk != "B.NAIVE") %>%
  left_join(source, by=c("COMBAT_ID_Time")) %>%
  ggplot(aes(Source_abrev, low_mut))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, lwd=0) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.5)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  facet_wrap(~pseudobulk, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  xlab("")+
  ylab("Proportion of repertoire <4 mutations")


#Median SMH replacement per pseudobulk ----
rsmut <- data.table::setDT(Heavy)[,list(Mean.r.mut=mean(r_mut), Max.r.mut=max(r_mut), Min.r.mut=min(r_mut), Median.r.mut=as.numeric(median(r_mut)), Std.r.mut=sd(r_mut),Mean.s.mut=mean(s_mut), Max.s.mut=max(s_mut), Min.s.mut=min(s_mut), Median.s.mut=as.numeric(median(s_mut)), Std.s.mut=sd(s_mut)), by=c("COMBAT_ID_Time", "pseudobulk")]

rsmut %>%
  dplyr::left_join(source, by=c("COMBAT_ID_Time")) %>%
  ggplot(aes(Median.s.mut, Median.r.mut, color=Source_abrev))+
  geom_point( size=0.8)+
  facet_grid(vars(pseudobulk), vars(Source_abrev))+
  geom_smooth(method="lm")+
  scale_color_manual(values=group.colors)

#Median rs mut
rsmut <- data.table::setDT(Heavy)[,list(Mean.r.mut=mean(r_mut), Max.r.mut=max(r_mut), Min.r.mut=min(r_mut), Median.r.mut=as.numeric(median(r_mut)), Std.r.mut=sd(r_mut),Mean.s.mut=mean(s_mut), Max.s.mut=max(s_mut), Min.s.mut=min(s_mut), Median.s.mut=as.numeric(median(s_mut)), Std.s.mut=sd(s_mut)), by=c("COMBAT_ID_Time", "pseudobulk")]

##junction length----
##bootstrap junction length function
getjl <- function(input, repertoire_id, cell_subset, B = 1000, subsample = 5, conf.level = 0.95) {
  data <- as.numeric(input$junction_length_HC)
  ID <- input[,c(repertoire_id, cell_subset)] #this is just to set the IDs later. If you add a new variable to your group split above, you will need to add the same variable at the end of here e.g. if I add cluster_id then that needs to go after "Source"
  options(`mosaic:parallelMessage` = FALSE)
  gini.boot <- mosaic::do(B) * mosaic::resample(data, subsample=5, replace = FALSE)
  names(gini.boot) <- "result"
  #As Gini can only accept a vector of more than two, when we only detect one clone, it creates a NaN by coercian. 
  #To fix this, we will manually reset NaNs to 1
  boot.mean <- mosaic::mean(~result, data = gini.boot)
  boot.se <- mosaic::sd(~result, data = gini.boot)
  boot.median <- mosaic::median(~result, data = gini.boot)
  g <- as.numeric(1 - conf.level)
  boot.lci <- mosaic::qdata(~result, g/2,data = gini.boot)
  boot.uci <- mosaic::qdata(~result,1-g, data = gini.boot)
  boot.stats <- data.table::data.table(mean.junction=boot.mean, se=boot.se, median.junction=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
  boot.stats <- cbind(ID, boot.stats)
  return(boot.stats)
}

cdrl_list <- sampletest %>% 
  select(COMBAT_ID_Time, Source_abrev, barcode_id, junction_length_HC, pseudobulk) %>%
  group_by(COMBAT_ID_Time, pseudobulk) %>%
  group_split()


cdrl.boot <-  mclapply(cdrl_list, function(x){
  getjl(input=x, repertoire_id = "COMBAT_ID_Time", cell_subset = "pseudobulk", subsample=5)#write subsample as a parameter
},mc.cores=5)

cdrl.boot <- data.table::rbindlist(cdrl.boot)

readr::write_csv(cdrl.boot, "out/cdrl.boot.csv")

cdrl.boot <- cdrl.boot %>%
  left_join(Heavy %>% select(COMBAT_ID_Time, Source_abrev) %>% unique, by=c("COMBAT_ID_Time"))
cdrl.boot$Source_abrev <- factor(cdrl.boot$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm", "Sepsis"))

cdrl.boot.dunn <- cdrl.boot %>%
  group_by(COMBAT_ID_Time, pseudobulk, Source_abrev) %>% 
  summarise(mean = mean(mean.junction, na.rm = TRUE)) %>%
  ungroup %>%
  group_by(pseudobulk) %>%
  rstatix::dunn_test(mean~Source_abrev, p.adjust.method = "fdr") %>%
  rstatix::add_xy_position(x="Source_abrev")
  
#plot
cdrl.boot.p <- cdrl.boot%>%
  group_by(COMBAT_ID_Time, pseudobulk, Source_abrev) %>% 
  summarise(mean = mean(mean.junction, na.rm = TRUE)) %>%
  ungroup %>%
  group_by(Source_abrev, pseudobulk) %>%
  summarise(mean.1 = mean(mean, na.rm = TRUE), 
            sd = sd(mean, na.rm = TRUE),
  n = n(), 
  lower.ci = quantile(mean, 0.25),
  upper.ci = quantile(mean, 0.75)) 
         

cdrl.boot.p$pseudobulk <- factor(cdrl.boot.p$pseudobulk, levels = c("B.NAIVE", "B.INT", "B.MEM", "PB"))
cdrl_p1<- cdrl.boot.p %>%
  ggplot(aes(pseudobulk, mean.1, color=Source_abrev))+
  geom_ribbon(aes(ymin=lower.ci, ymax=upper.ci, group=Source_abrev, fill=Source_abrev),alpha=0.2, color=NA)+
  geom_line(aes(group=Source_abrev))+
  scale_color_manual(values=group.colors)+
  scale_fill_manual(values=group.colors)+
  facet_wrap(~Source_abrev)+
  xlab("")+
  ylab("Mean junction length (boot=1000)")

pdf("COMBAT_figure_plots/cdrl_plot.pdf", height=5, width=7)
cdrl_p1
dev.off()

cdrl.plot <- cdrl.boot %>%
  filter(pseudobulk=="PB")%>%
  mutate(anno= "Plasmablast")%>%
  group_by(COMBAT_ID_Time, pseudobulk, Source_abrev) %>% 
  summarise(mean = mean(mean.junction, na.rm = TRUE)) %>%
  ggplot(aes(Source_abrev, mean))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, lwd=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=0.5)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  facet_wrap(~pseudobulk, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  xlab("")+
  ylab("Mean junction length (boot=1000)")+
  theme(legend.position = "none")

cdrl.boot.dunn %>%filter(p.adj.signif != "ns")
pdf("COMBAT_figure_plots/cdrl_pb.pdf", height=4, width=2)
cdrl.plot + ggpubr::stat_pvalue_manual(cdrl.boot.dunn %>% filter(pseudobulk=="PB")%>%filter(p.adj.signif != "ns")  , label = "p.adj.signif", tip.length = 0.01, y.position = 22, step.increase = 0.07)
dev.off()

##Bulk repertoire IGHV4-34 analysis ----
raw <- data.table::fread("/Users/bosun/Documents/GitHub/COMBAT_rep/data/All_raw_values_COMBAT_BCR1_FINAL.txt") %>% filter(!grepl("Batch", Source))
colnames(raw)
raw$RNASeq_sample_ID %in% Heavy$RNASeq_sample_ID %>% summary
## stats for unique BCRs for igha1 2 and ighg1
iso <- raw[,c(2,3,10,52:60)]

iso_list <- iso %>%
  reshape2::melt(id.vars=c("RNASeq_sample_ID", "Source", "COMBAT_ID")) %>%
  filter(value >= 0)%>%
  filter(!grepl("IGHDM", variable)) %>%
  filter(!grepl("class_switched", variable)) %>%
  filter(!grepl("IGHD", variable)) %>%
  left_join(source2, by=c("Source"))

iso_list$Source_abrev <- factor(iso_list$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))

iso_list <- split(iso_list, as.character(iso_list$variable))
names(iso_list) <- c("Unmutated V4-34 AVY/NHS IGHA", "Unmutated V4-34 AVY/NHS IGHA2", "Unmutated V4-34 AVY/NHS IGHG1","Unmutated V4-34 AVY/NHS IGHG2", "Unmutated V4-34 AVY/NHS IGHG3", "Unmutated V4-34 AVY/NHS IGHM")
iso_names <- names(iso_list)
iso_labels <- as.list(c("Unmutated V4-34 AVY/NHS IGHA", "Unmutated V4-34 AVY/NHS IGHA2", "Unmutated V4-34 AVY/NHS IGHG1","Unmutated V4-34 AVY/NHS IGHG2", "Unmutated V4-34 AVY/NHS IGHG3", "Unmutated V4-34 AVY/NHS IGHM"))
names(iso_labels) <- names(iso_list)

iso_list <- lapply(iso_names, function(x){
  y <- iso_list[[x]] %>%
    mutate(variable=iso_labels[[x]])
  y$Source_abrev <- factor(y$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))
return(y)
})

names(iso_list) <- c("Unmutated V4-34 AVY/NHS IGHA", "Unmutated V4-34 AVY/NHS IGHA2", "Unmutated V4-34 AVY/NHS IGHG1","Unmutated V4-34 AVY/NHS IGHG2", "Unmutated V4-34 AVY/NHS IGHG3", "Unmutated V4-34 AVY/NHS IGHM")

iso_dunn <- lapply(iso_names, function(x) {
  iso_list[[x]] %>%
    rstatix::dunn_test(value ~ Source_abrev, p.adjust.method = "fdr")%>%
  add_x_position(x="Source") %>%
    mutate(variable=iso_labels[[x]])
})

names(iso_dunn) <- names(iso_list)

iso_plotlist <- lapply(iso_names, function(x){
  iso_boxplot <- iso_list[[x]] %>%
    ggplot(aes(Source_abrev, value))+
    geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
    geom_jitter(aes(color=Source_abrev),width=0.1, size=0.6)+
    stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.6, alpha=90)+
    #scale_fill_manual(values=group.colors)+
    #scale_color_manual(values=group.colors)+
    xlab("")+
    ylab("Percentage of isotype repertoire")+
    theme(legend.position = "none")+
    facet_wrap(~variable)+
    scale_color_manual(values=group.colors)+
    scale_fill_manual(values=group.colors)+
    theme(plot.title = element_text(size = 7, face = "bold"))
  plot<- iso_boxplot + ggpubr::stat_pvalue_manual(iso_dunn[[x]] %>%
                                                    filter(p.adj.signif !="ns")  , label = "p.adj.signif", tip.length = 0.01, y.position= (max(iso_list[[x]]$value) + max(iso_list[[x]]$value)*0.1),step.increase = 0.06, hide.ns = TRUE)
  
return(plot)  
})

patchwork::wrap_plots(iso_plotlist)

pdf("COMBAT_figure_plots/434_avy_nhs_plot_2.pdf", height= 8, width=8)
patchwork::wrap_plots(iso_plotlist)
dev.off()

##Bulk repertoire class switch plotting NEED TO FILL----
##looking at clonal binning----
#for clonal binning, rather than look at mutations which are not normalised for length. I will bin by V identity
Heavy <- Heavy %>% 
  mutate(low_mut = ifelse(v_identity >=95, "low_mut", "high_mut")) %>%
  mutate(low_exp = ifelse(v_identity >=95 & clonal_abundance >=2,">=95% IGHV identity expanded", 
                          ifelse(v_identity >=95 & clonal_abundance <2, ">=95% IGHV identity unexpanded", 
                                 ifelse(v_identity <95 & clonal_abundance >=2, "<95% IGHV identity expanded", 
                                        ifelse(v_identity <95 & clonal_abundance <2, "<95% IGHV identity unexpanded", "not_assigned")))))


nums <-  Heavy %>%
  select(COMBAT_ID_Time,Source_abrev, low_exp) %>%
  unique() %>%
  group_by(Source_abrev,low_exp) %>%
  add_tally(name="number_of_individuals") %>%
  select(-COMBAT_ID_Time) %>%
  unique() %>%
  as.data.frame()

prop_bins_per_pseudobulk <- Heavy %>% 
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  group_by(COMBAT_ID_Time, pseudobulk, low_exp) %>% 
  add_tally(name="sum_clonal") %>%
  dplyr::select(COMBAT_ID_Time, Source_abrev, pseudobulk, low_exp,sum,sum_clonal) %>% 
  unique %>%
  mutate(prop=sum_clonal/sum) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  dplyr::select(-sum, -sum_clonal) %>%
  ungroup %>%
  unique %>%
  spread(low_exp, value=prop)

prop_bins_per_pseudobulk[is.na(prop_bins_per_pseudobulk)] <- 0

prop_bins_per_pseudobulk <- prop_bins_per_pseudobulk %>%
  reshape2::melt(id.vars=c("COMBAT_ID_Time", "Source_abrev", "pseudobulk"))

prop_bins_per_pseudobulk$Source_abrev <- factor(prop_bins_per_pseudobulk$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Sepsis"))


p1 <- prop_bins_per_pseudobulk %>%
  ggplot(aes(Source_abrev, value, color=Source_abrev))+
  geom_boxplot(color="black", outlier.shape = NA)+
  geom_jitter(width=0.2, size=0.5)+
  facet_grid(vars(pseudobulk), vars(variable))+
  scale_color_manual(values=group.colors)+
  ylab("Proportion of clonal bins per pseudobulk")


prop_bins_per_pseudobulk_list<- prop_bins_per_pseudobulk%>%
  filter(grepl("PB", pseudobulk)) %>%
  filter(!(grepl("Flu", Source_abrev))) %>%
  group_by(variable) %>%
  group_split()

names(prop_bins_per_pseudobulk_list) <- unique(prop_bins_per_pseudobulk$variable)


z <- lapply(prop_bins_per_pseudobulk_list, function(y){
    y.list <- split(y, y$pseudobulk)
  
   dunnlist <-lapply(y.list, function(z){
     if(sum(z$value) != 0){
     dunn<- rstatix::dunn_test(z, value ~ Source_abrev, p.adjust.method = "fdr",detailed = FALSE) %>%
       rstatix::add_xy_position(x="Source_abrev") %>%
       mutate(pseudobulk=unique(z$pseudobulk)) %>%
       mutate(variable=unique(z$variable))
     return(dunn)  
     }else{
       NULL
     }
   })
    return(dunnlist)  
})


p1

zbind <- lapply(z, function(x){
  x[sapply(x, is.null)] <- NULL
  data.table::rbindlist(x) 
})

binlist <- names(prop_bins_per_pseudobulk_list)

plotlist_bins<- lapply(binlist, function(x){
  prop_list <- prop_bins_per_pseudobulk_list[[x]] 
  
  pbprop<- prop_list %>%
    ggplot(aes(Source_abrev, value, color=Source_abrev))+
    geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
    geom_jitter(aes(color=Source_abrev),width=0.1, size=0.6)+
    stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.6, alpha=90)+
    scale_fill_manual(values=group.colors)+
    scale_color_manual(values=group.colors)+
    facet_wrap(~variable)+
    ylab("Proportion clonal bins within Plasmablasts")+
    theme(legend.position="none")+
    xlab("")
  
  if(nrow(zbind[[x]] %>% filter(p.adj.signif != "ns"))!=0){
  p1 <- pbprop+ggpubr::stat_pvalue_manual(zbind[[x]] %>% filter(p.adj.signif != "ns") , label = "p.adj.signif", tip.length = 0.01, y.position= (max(prop_list$value) + max(prop_list$value)*0.1), step.increase = 0.1,)
  
  return(p1)
  } else {
    return(pbprop)
  }
})

pbprob <- patchwork::wrap_plots(plotlist_bins, nrow=2)
pdf("COMBAT_figure_plots/clonal_bins_pb.pdf", height=8, width=5.3)
pbprob
dev.off()

##clonal overlaps per pseudobulk
Heavy$pseudobulk <- factor(Heavy$pseudobulk, levels=c("B.NAIVE", "B.INT","B.MEM",  "PB"))

sourcesplit <- split(Heavy, Heavy$Source_abrev)

overlaplist <- lapply(sourcesplit, function(x){
  y<- table(x$clone_per_replicate, x$pseudobulk)
  y[(y>1)] <- 1
  
    mat <- crossprod(as.matrix(y))            # calculate the overlap               
    #mat <- floor((mat * 100 / diag(mat))) 
    mat[upper.tri(mat)] <- NA
    mat[(mat==0)] <- NA
    
  return(mat)
})

overlapnames <- names(overlaplist)
overlaplong <- lapply(overlapnames, function(x){
  z <- overlaplist[[x]]%>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    dplyr::rename(`Number of clonal overlaps`=value) %>%
    mutate(Source= paste(x))
})

overlap <- data.table::rbindlist(overlaplong)
overlap$Source <- factor(overlap$Source, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))

clonal_overlap_pseudo <- overlap %>%
  ggplot(aes(Var1, Var2, fill=`Number of clonal overlaps`))+
  geom_tile()+
  geom_text(aes(label=`Number of clonal overlaps`))+
  scale_fill_distiller(palette = "Spectral", , trans = "log10")+
  #viridis::scale_fill_viridis()+
  theme(legend.position ="right")+
  facet_wrap(~Source)+
  xlab("")+
  ylab("")

pdf("COMBAT_figure_plots/clonal_overlap_pseudo.pdf", height=6, width=9)
clonal_overlap_pseudo
dev.off()

##class and clonal overlap supplementary----
c_genetest <- Heavy
c_genetest$c_gene <- factor(c_genetest$c_gene , levels= c("IGHM", "IGHD", "IGHG3","IGHG1", "IGHA1","IGHG2","IGHG4","IGHE","IGHA2", "None"))

sourcesplit <- split(c_genetest, c_genetest$Source_abrev)

c_overlaplist <- lapply(sourcesplit, function(x){
  y<- table(x$clone_per_replicate, x$c_gene)
  y[(y>1)] <- 1
  
  mat <- crossprod(as.matrix(y))            # calculate the overlap               
  #mat <- floor((mat * 100 / diag(mat))) 
  mat[upper.tri(mat)] <- NA
  mat[(mat==0)] <- NA
  
  return(mat)
})

c_overlapnames <- names(c_overlaplist)
c_overlaplong <- lapply(c_overlapnames, function(x){
  z <- c_overlaplist[[x]]%>%
    reshape2::melt() %>%
    filter(!is.na(value)) %>%
    dplyr::rename(`Number of clonal overlaps`=value) %>%
    mutate(Source= paste(x))
})

c_overlap <- data.table::rbindlist(c_overlaplong)
c_overlap$Source <- factor(c_overlap$Source, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))

clonal_overlap_c_gene<- c_overlap %>%
  filter(!grepl("IGHE|None", Var1)) %>%
  filter(!grepl("IGHE|None", Var2)) %>%
  ggplot(aes(Var1, Var2, fill=`Number of clonal overlaps`))+
  geom_tile()+
  geom_text(aes(label=`Number of clonal overlaps`))+
  scale_fill_distiller(palette = "Spectral", , trans = "log10")+
  #viridis::scale_fill_viridis()+
  theme(legend.position ="right")+
  facet_wrap(~Source)+
  xlab("")+
  ylab("")

pdf("COMBAT_figure_plots/clonal_overlap_c_gene.pdf", height=6, width=12)
clonal_overlap_c_gene
dev.off()

##Plot convergent clones per group - not selected for those just convergent in COVID----
convergent_clones <- Heavy %>% 
  select(Source_abrev, clone_global) %>%
  unique %>%
  group_by(clone_global) %>%
  add_tally() %>%
  filter(n>1) %>%
  arrange(clone_global) %>%
  mutate(Source_abrev = gsub("CS|CM|CC|CComm", "COVID19", Source_abrev))

convergent_clones <- convergent_clones %>%
  group_by(clone_global) %>%
  summarise_all(funs(paste0(unique(.[!is.na(.)]), collapse= "_"))) %>%
  dplyr::rename(overlap_labels= Source_abrev) %>%
  mutate(conv_label = "convergent")

convergent <- Heavy %>%
  select(COMBAT_ID_Time, barcode_id, pseudobulk, junction_aa, v_gene_HC, j_gene_HC, Source_abrev, clone_global, umis) %>%
  left_join(convergent_clones, by=c("clone_global")) %>%
  mutate(conv_label = ifelse(is.na(conv_label), "non_convergent", conv_label))%>%
  mutate(overlap_labels = ifelse(is.na(overlap_labels), "non_convergent", overlap_labels))

convergent_prop <- convergent %>%
  group_by(COMBAT_ID_Time, pseudobulk)%>%
  add_tally(name = "sum") %>%
  group_by(COMBAT_ID_Time, pseudobulk, overlap_labels) %>% 
  add_tally(name="sum_p") %>%
  select(COMBAT_ID_Time, Source_abrev, pseudobulk,overlap_labels,sum,sum_p) %>% 
  unique %>%
  mutate(conv_prop=sum_p/sum) %>%
  arrange(pseudobulk,COMBAT_ID_Time) %>%
  select(-sum, -sum_p) %>%
  spread(overlap_labels, value=conv_prop)

convergent_prop[is.na(convergent_prop)] <- 0

convergent_prop <- convergent_prop %>%
  reshape2::melt(id.vars=c("COMBAT_ID_Time", "Source_abrev", "pseudobulk"))
convergent_prop$Source_abrev <- factor(convergent_prop$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))

convergent_p1 <- convergent_prop %>%
  filter(variable != "non_convergent") %>%
  mutate(variable = gsub("HV_COVID19","COVID19_HV", variable)) %>%
  mutate(variable = gsub("Sepsis_COVID19", "COVID19_Sepsis", variable)) %>%
  #mutate(overlap_labels=gsub("CS|CM|CC|CComm", "COVID19", overlap_labels)) %>%
  #mutate(overlap_labels=gsub("COVID19_COVID19", "COVID19", overlap_labels)) %>%
  ggplot(aes(Source_abrev, value+0.01, fill=Source_abrev))+
  geom_violin(scale="width", alpha=0.8)+
  geom_boxplot(color="black", fill="white", width=0.1, alpha=0.5, outlier.shape = NA)+
  geom_jitter(width=0.2, alpha=0.5)+
  facet_grid(vars(variable), vars(pseudobulk))+
  scale_fill_manual(values=group.colors)+
  scale_y_log10()

pdf("COMBAT_figure_plots/convergent_pseudobulks.pdf", height=5, width=7)
convergent_p1
dev.off()

##Bulk BCR class switching plot ----
file = "data/Class-switching COMBAT/Max_severity_statistical_summary_age adjusted.txt"
p <- as.matrix(read.csv(file, head=TRUE, sep="\t"))
p=p[which(p[,"type"]=="Relative_class_switching_absolute"),]
headers = colnames(p)
classes = gsub("Relative_class_switching_absolute..","",p[,"class1"])
classes = gsub("IGHD.M","IGHD/M", classes,fixed = T)

p_head = headers[grep("pvalues.", headers)]
m_head = headers[grep("means.", headers)]
p_values = matrix(data = -2, nrow = length(classes), ncol = length(p_head), dimnames = c(list(classes),list(p_head)))
means = matrix(data = -2, nrow = length(classes), ncol = length(m_head), dimnames = c(list(classes),list(m_head)))
for(i in c(1:length(m_head))){means[, m_head[i]] = as.numeric(p[, m_head[i]])}
for(i in c(1:length(p_head))){p_values[, p_head[i]] = as.numeric(p[, p_head[i]])}

class = strsplit(classes,".",fixed = T)
class1 = NULL
class2 = NULL
for(i in c(1:length(class))){
  class1 = c(class1, class[[i]][1])
  class2 = c(class2, class[[i]][2])
}

match_edges = rep(-1, length(class1))
for(i in c(1:length(class))){
  w = intersect(which(edge_ids[,1]== class1[i]),which(edge_ids[,2]== class2[i]))
  if(length(w)==0){w = intersect(which(edge_ids[,1]== class2[i]),which(edge_ids[,2]== class1[i]))}
  if(length(w)==1){match_edges[i]=w}}


df <- as.data.frame(p) %>% 
  mutate(From= c("IGHA1", "IGHD/M", "IGHG1", "IGHA1", "IGHG3", "IGHD/M", "IGHG1", "IGHG2", "IGHD/M", "IGHD/M", "IGHD/M", "IGHG1", "IGHG3"))%>% 
  mutate(type=From) %>%
  mutate(To= c("IGHA2", "IGHA1", "IGHA1", "IGHG2", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG1", "IGHG2", "IGHG3", "IGHG2", "IGHG2"))

df <- df %>% select(From, To, pvalues.COVID_MILD, pvalues.COVID_SEV, 
                    pvalues.COVID_CRIT, pvalues.COVID_HCW_MILD, 
                    pvalues.Sepsis, means.HV,means.COVID_MILD, 
                    means.COVID_SEV, means.COVID_CRIT, 
                    means.COVID_HCW_MILD, means.Sepsis, p_value_HV_COV_all) %>%
  mutate(From_to=paste(From, To, sep="_"))
df <- df %>% mutate(pvalues.HV=p_value_HV_COV_all)

group_id <- c("HV", "COVID_MILD", "COVID_SEV","COVID_CRIT", "COVID_HCW_MILD", "Sepsis")

el.g <- get.edgelist(g) %>% as.data.frame() %>% mutate(From_to=paste(V1, V2, sep="_"))
el.main <- get.edgelist(g) %>% as.data.frame() 
el.main.g <- igraph::graph_from_edgelist(el.main %>% as.matrix())
el.main.n <- asNetwork(el.main.g)
el.main.n1 <- ggnetwork(el.main.n, layout = "circle", arrow.gap=0.12)
el.main.n1$vertex.names <- factor(el.main.n1$vertex.names, levels = c("IGHD/M", "IGHG3","IGHG1", "IGHA1","IGHG2","IGHG4","IGHE","IGHA2"))
plot_n1 <- ggplot(el.main.n1, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend), color = "white",arrow = arrow(length = unit(6, "pt"), type="closed")) +
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend),arrow = arrow(length = unit(6, "pt"), type="closed"), linejoin='mitre') +
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend), color="white", linejoin='mitre') +
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend),linetype="11", linejoin='mitre') +
  scale_color_manual(values=c("gray", "red"))+
  new_scale_color()+
  geom_nodes(aes(color = vertex.names), size=22, alpha=0.3) +
  geom_nodes(size=18, color="white") +
  geom_nodes(aes(color = vertex.names), size=18, alpha=0.5) +
  ggsci::scale_color_d3()+
  geom_text(aes(label=vertex.names))+ 
  theme_blank()+
  xlim(-0.2,1.2)+
  ylim(-0.1, 1.1)+ 
  scale_size(range = c(1, 2))+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

plot_list <- list()
for (i in 1: length(group_id)){
  subset2 <- el.g %>% left_join(df, by="From_to")
  subset2[,11:16][is.na(subset2[,11:16])] <- 0
  df2 <- subset2 %>% select(V1,V2) %>% as.matrix() 
  f1 <- igraph::graph_from_edgelist(df2)
  E(f1)$weight <- as.numeric(subset2[[paste("means",group_id[[i]], sep=".")]])
  E(f1)$p_sig <- subset2[[paste("pvalues",group_id[[i]], sep=".")]]
  f2 <- asNetwork(f1)
  n1 <- ggnetwork(f2, layout = "circle", arrow.gap=0.15)
  n1$weight <- as.numeric(n1$weight)
  n1$p_sig <- as.numeric(n1$p_sig)
  n1 <- n1 %>% arrange(desc(p_sig))
  n1$vertex.names <- factor(n1$vertex.names, levels = c("IGHD/M", "IGHG3","IGHG1", "IGHA1","IGHG2","IGHG4","IGHE","IGHA2"))
  n1 <- n1 %>% mutate(weight= ifelse(weight!= 0, weight, NA))
  if(group_id[[i]]=="Sepsis"|group_id[[i]]=="HV"){
    p1 <- ggplot(n1, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight), color = "white",arrow = arrow(length = unit(6, "pt"), type="closed")) +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight, color=p_sig<0.05, alpha=-p_sig),arrow = arrow(length = unit(6, "pt"), type="closed"), linejoin='mitre') +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight), color="white", linejoin='mitre') +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight, color=p_sig<0.05, alpha=-p_sig),linetype="11", linejoin='mitre') +
      scale_color_manual(values=c("gray", "red"))+
      new_scale_color()+
      geom_nodes(aes(color = vertex.names), size=22, alpha=0.3) +
      geom_nodes(size=18, color="white") +
      geom_nodes(aes(color = vertex.names), size=18, alpha=0.5) +
      ggsci::scale_color_d3()+
      geom_text(aes(label=vertex.names))+ 
      theme_blank()+
      xlim(-0.2,1.2)+
      ylim(-0.1, 1.1)+ 
      scale_size(range = c(1, 2))+
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
      ggtitle(paste0(group_id[[i]]))
  }else{
    p1 <- ggplot(n1, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight), color = "white",arrow = arrow(length = unit(6, "pt"), type="closed")) +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight, color=p_sig<0.05, alpha=-p_sig),arrow = arrow(length = unit(6, "pt"), type="closed"), linejoin='mitre') +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight), color="white", linejoin='mitre') +
      geom_edges(aes(x = x, y = y, xend = xend, yend = yend, size=weight, color=p_sig<0.05, linetype=p_sig>0.05, alpha=-p_sig), linejoin='mitre') +
      scale_color_manual(values=c("gray", "red"))+
      new_scale_color()+
      geom_nodes(aes(color = vertex.names), size=22, alpha=0.3) +
      geom_nodes(size=18, color="white") +
      geom_nodes(aes(color = vertex.names), size=18, alpha=0.5) +
      ggsci::scale_color_d3()+
      geom_text(aes(label=vertex.names))+ 
      theme_blank()+
      xlim(-0.2,1.2)+
      ylim(-0.1, 1.1)+ 
      scale_size(range = c(1, 2))+
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5))+
      ggtitle(paste0(group_id[[i]]))
  }
  plot_list[[i]] <- p1
}

pdf("COMBAT_figure_plots/bulk_class_switch_network.pdf", height= 4, width=25)
patchwork::wrap_plots(plot_list, nrow=1)
dev.off()


##COMBAT clonal overlap upset plot----
#This captures shared clones by similarity and clonal relatedness not just 100% identity
shared <- Heavy %>% select(clone_global, Source_abrev, pseudobulk, barcode_id, COMBAT_ID_Time) %>% filter(Source_abrev !="Flu")%>% group_by(clone_global) %>% add_tally() %>% arrange(clone_global)
VENN.LIST <- shared %>% select(Source_abrev, pseudobulk,clone_global) %>% split(., .$Source_abrev, drop = TRUE)
VENN.LIST <- lapply(VENN.LIST, function(x){x <- x$clone_global})
upsetdf<- UpSetR::fromList(VENN.LIST)
a.no0 = upsetdf[ rowSums(upsetdf)!=1, ] 
names(VENN.LIST)

pdf("COMBAT_figure_plots/clonal_upsetplot.pdf", height= 5, width=10)
UpSetR::upset(a.no0, sets = rev(c("HV", "CM","CS", "CC","CComm", "Sepsis")) , mb.ratio = c(0.55, 0.45), order.by = "freq", 
              keep.order = TRUE)
dev.off()


##Sequence logos for COMBAT overlapping clones----
nonuniqueconv_names<- Heavy %>%
  filter(Source_abrev != "Flu")%>%
  select(barcode_id, Source_abrev, clone_global, junction_aa, total_mut) %>%
  select(Source_abrev, clone_global) %>%
  unique %>%
  group_by(clone_global) %>%
  add_tally() %>%
  select(Source_abrev, clone_global, n) %>%
  filter(n>1) %>%
  arrange(desc(n)) %>%
  filter(n==4)

conv_names<- Heavy %>%
  filter(Source_abrev != "Flu")%>%
  filter(Source_abrev != "HV")%>%
  filter(Source_abrev != "Sepsis")%>%
  select(barcode_id, Source_abrev, clone_global, junction_aa, total_mut) %>%
  select(Source_abrev, clone_global) %>%
  unique %>%
  group_by(clone_global) %>%
  add_tally() %>%
  select(Source_abrev, clone_global, n) %>%
  filter(n>1) %>%
  arrange(desc(n))

conv_clones <- Heavy[c(Heavy$clone_global %in% conv_names$clone_global),] %>% select(barcode_id,Source_abrev, clone_global, junction_aa) %>% filter(Source_abrev != "Flu")
conv_clones %>%
  select(clone_global) %>% unique

conv_clones <- split(conv_clones, conv_clones$clone_global)

library(Biostrings)
aalist <- lapply(conv_clones, function(x){
  y <- AAStringSet(x$junction_aa)
  names(y) <- paste(x$clone_global, x$Source_abrev, sep="_")
  return(y)
})

msalist <- lapply(aalist, function(x){
  x <- msa(x)
  return(x)
})
length(unmasked(msalist[[8]]))



library(ggseqlogo)

seqlist <- lapply(msalist, function(x){
  y <- as.vector(unmasked(x))
  return(y)
})

seqlist_l <- lapply(seqlist, function(x){
  x <- length(x)
})

unlist(seqlist_l) %>%as.data.frame %>% arrange(desc(.))

inter_4 <- ggseqlogo::ggseqlogo(seqlist$`1266-ALLCOMBAT`)
inter_3 <- ggseqlogo::ggseqlogo(seqlist$`25767-ALLCOMBAT`)

pdf("COMBAT_figure_plots/intersecting_cluster_2.pdf", height = 1.7, width = 3.5)
ggseqlogo::ggseqlogo(seqlist$`1266-ALLCOMBAT`)
dev.off()
pdf("COMBAT_figure_plots/intersecting_cluster_1.pdf", height = 1.7, width = 3.5)
ggseqlogo::ggseqlogo(seqlist$`25767-ALLCOMBAT`)
dev.off()



#get networks for all clones----
dist_list <- lapply(conv_clones, function(x){
  dist <- stringdistmatrix(x$junction_aa, x$junction_aa, method = "lv") 
  rownames(dist) <- x$barcode_id
  colnames(dist) <- x$barcode_id
  dist[upper.tri(dist)] <- NA
  dist <- reshape2::melt(dist)
  dist <- dist %>% mutate(value = value+1) %>% filter(!is.na(value))
  
  return(dist)
})

dist_bound <- data.table::rbindlist(dist_list)

g <- igraph::graph_from_edgelist(dist_bound %>% select(Var1, Var2) %>% as.matrix(), directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
n <- asNetwork(g)
g <- ggnetwork(n)

g <- g %>% unique
g <- g %>% left_join(Heavy %>% select(barcode_id, c_gene, pseudobulk, Source_abrev, clone_global), by=c("vertex.names" ="barcode_id"))


ggplot(g, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color="black", size=0.05) +
  geom_nodes(aes(color=Source_abrev), size=0.5)+
  theme_blank()+
  scale_color_manual(values=group.colors)+
  guides(linetype=FALSE,
         colour = guide_legend(override.aes = list(size=2))) #Remove legends

##Pan BCRs for published SARS-CoV2-specific mAbs (http://opig.stats.ox.ac.uk/webapps/covabdab/)----
covdb <- data.table::fread("data/CoV-AbDab_120121.csv")

#filter for human seqs
hucovdb <- covdb %>%
  filter(grepl("Human|human", Origin)) %>%
  filter(grepl("Human|human", `Heavy V Gene`)) 

hucovdb <- hucovdb %>% 
  mutate(v_gene_hcov = gsub("\\ .*", "", `Heavy V Gene`)) %>%
  mutate(j_gene_hcov = gsub("\\ .*", "", `Heavy J Gene`)) %>%
  mutate(v_gene_lcov = gsub("\\ .*", "", `Light V Gene`)) %>%
  mutate(j_gene_lcov = gsub("\\ .*", "", `Light J Gene`)) %>%
  filter(CDRH3 != "ND")

Heavy$hcdr3 <- gsub("^C|W$","", Heavy$junction_aa)

covdist <- stringdistmatrix(hucovdb$CDRH3,Heavy$hcdr3, method="lv")
rownames(covdist) <- hucovdb$Name
colnames(covdist) <- Heavy$barcode_id
covdist[upper.tri(covdist)] <- NA
covdistm <- reshape2::melt(covdist)
colnames(covdistm) <- c("seq_from", "seq_to", "lv_dist")

metacovdist <- covdistm %>%
  filter(lv_dist <=2)%>%
  left_join(hucovdb %>% select(Name, `Neutralising Vs`,`Protein + Epitope`,v_gene_hcov,j_gene_hcov, CDRH3, v_gene_lcov, j_gene_lcov, CDRL3), by=c("seq_from"="Name")) %>%
  left_join(Heavy %>% select(barcode_id,v_gene_HC, j_gene_HC, hcdr3, COMBAT_ID_Time, Source_abrev, pseudobulk, clone_global), by=c("seq_to"="barcode_id")) %>%
  left_join(Light %>% select(barcode_id,v_gene_LC, j_gene_LC, junction_aa), by=c("seq_to"="barcode_id"))

metacovdist[c(metacovdist$v_gene_hcov == metacovdist$v_gene_HC),] 

hucovdb$concat <- paste(hucovdb$v_gene_hcov, hucovdb$j_gene_hcov, sep="|")
hucovdb$j_length <- nchar(hucovdb$CDRH3)


vjmatch <- Heavy[c(Heavy$VJ_HC %in% hucovdb$concat),] 

vjlength <- vmatch %>%
  left_join(hucovdb %>% select(Name, concat, CDRH3, j_length), by=c("VJ_HC"="concat")) %>% 
  filter(nchar(hcdr3) == j_length) 

covdist <- stringdistmatrix(vjlength$CDRH3,vjlength$hcdr3, method="lv")
rownames(covdist) <- vjlength$Name
colnames(covdist) <- vjlength$barcode_id
covdist[upper.tri(covdist)] <- NA
covdistm <- reshape2::melt(covdist)
colnames(covdistm) <- c("seq_from", "seq_to", "lv_dist")

covmatch <- covdistm %>%
  filter(lv_dist <=3) %>%
  left_join(hucovdb %>% select(Name, `Neutralising Vs`,`Protein + Epitope`,v_gene_hcov,j_gene_hcov, CDRH3, v_gene_lcov, j_gene_lcov, CDRL3, Sources), by=c("seq_from"="Name")) %>%
  left_join(Heavy %>% select(barcode_id,v_gene_HC, j_gene_HC, hcdr3, COMBAT_ID_Time, Source_abrev, pseudobulk, clone_global), by=c("seq_to"="barcode_id")) %>%
  left_join(Light %>% select(barcode_id,v_gene_LC, j_gene_LC, junction_aa), by=c("seq_to"="barcode_id")) %>%
  select(seq_from, seq_to,lv_dist, Source_abrev, COMBAT_ID_Time, `Neutralising Vs`,`Protein + Epitope`, v_gene_HC,v_gene_hcov, j_gene_HC,j_gene_hcov, CDRH3, hcdr3, v_gene_LC, v_gene_lcov, j_gene_LC,j_gene_lcov, CDRL3, junction_aa, pseudobulk, Sources) %>%
  unique 

covmatch$junction_aa <- gsub("^C|F$","", covmatch$junction_aa)
covmatch <- covmatch %>% 
  filter(v_gene_HC == v_gene_hcov)

covmab_el.main.g <- igraph::graph_from_edgelist(covmatch %>% select(seq_from, seq_to) %>% as.matrix())
E(covmab_el.main.g)$weight <- as.numeric(covmatch$lv_dist)

covmab_el.main.n <- asNetwork(covmab_el.main.g)
covmab_el.main.n1 <- ggnetwork(covmab_el.main.n)
covmab_el.main.n1 <- covmab_el.main.n1 %>% left_join(covmatch %>% select(seq_from, `Protein + Epitope`), by=c("vertex.names"="seq_from"))
covmab_el.main.n1 <- covmab_el.main.n1 %>% left_join(covmatch %>% select(seq_to, Source_abrev), by=c("vertex.names"="seq_to"))
covmab_el.main.n1 <- covmab_el.main.n1 %>% left_join(covmatch %>% select(seq_from, Sources), by=c("vertex.names"="seq_from"))
covmab_el.main.n1$Sources <- gsub("\\,.*", ",", covmab_el.main.n1$Sources)
covmab_el.main.n1$Sources <- gsub("\\(.*", "", covmab_el.main.n1$Sources)

covmab_el.main.n1 <- covmab_el.main.n1 %>% mutate(`Protein + Epitope`=paste(`Protein + Epitope`, Sources, sep="  |  "))

covmab_el.main.n1 <- covmab_el.main.n1 %>% mutate(`Protein + Epitope`=ifelse(`Protein + Epitope`=="NA  |  NA", as.character(Source_abrev), `Protein + Epitope`) ) %>%
  select(-Source_abrev) 

covmab_el.main.n1 <- covmab_el.main.n1 %>% left_join(covmatch %>% select(seq_to, pseudobulk), by=c("vertex.names"="seq_to")) %>% unique
covmab_el.main.n1 <- covmab_el.main.n1 %>% dplyr::rename(lv_dist=weight)

cov_mab_p1 <- ggplot(covmab_el.main.n1, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(x = x, y = y, xend = xend, yend = yend, linetype=as.factor(lv_dist), color=lv_dist), linejoin='mitre') +
  scale_color_gradient(high="gray", low="red")+
  new_scale_color()+
  geom_nodes(data = covmab_el.main.n1 %>% filter(!is.na(pseudobulk)),aes(color = `Protein + Epitope`), size=15, alpha=0.4) +
  geom_nodes(color="white", size=9, alpha=1) +
  geom_nodes(data = covmab_el.main.n1 %>% filter(!is.na(pseudobulk)),aes(color = `Protein + Epitope`), size=9, alpha=0.7) +
  geom_nodes(data = covmab_el.main.n1 %>% filter(is.na(pseudobulk)),aes(color = `Protein + Epitope`), size=5, alpha=1) +
  geom_edgetext(aes(label = lv_dist), color = "grey25", size=2) +
  ggsci::scale_color_rickandmorty()+
  theme_blank()+
  xlim(-0.2,1.2)+
  ylim(-0.1, 1.1)+ 
  scale_size(range = c(1, 2))+
  geom_text(aes(label=pseudobulk), color="black", size=2, max.overlaps = Inf, 
            segment.size	=0.3)+
  guides(linetype=FALSE,
         colour = guide_legend(override.aes = list(size=3))) #Remove legends

pdf("COMBAT_figure_plots/cov_mab_network.pdf", height= 8, width=8)
cov_mab_p1
dev.off()

##Sequence logos for COMBAT/SARS-CoV2-specific mAb matches ----
covmatch_clones <- split(covmatch, covmatch$seq_to)

aalist <- lapply(covmatch_clones, function(x){
  y <- AAStringSet(cbind(x$CDRH3, x$hcdr3))
  names(y) <- paste(x$seq_to, x$Source_abrev, sep="_")
  return(y)
})

msalist <- lapply(aalist, function(x){
  x <- msa(x)
  return(x)
})

seqlist <- lapply(msalist, function(x){
  y <- as.vector(unmasked(x))
  return(y)
})

seqlist_l <- lapply(seqlist, function(x){
  x <- length(x)
})

logo_plots<- lapply(seqlist, function(x){
  ggseqlogo::ggseqlogo(x)
})

patchwork::wrap_plots(logo_plots)

lapply(names(logo_plots), function(x) ggsave(filename=paste("COMBAT_figure_plots/",x,"_covmAb_match",".pdf",sep=""), plot=logo_plots[[x]], height = 1.7, width = 3.5))




