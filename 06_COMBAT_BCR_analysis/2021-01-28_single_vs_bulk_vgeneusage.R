##Bulk BCR v gene analysis - supplementary
##Author: Bo Sun
##Group: Bashford-Rogers

library(MKmisc)
library(tidyverse)
library(ggplot2); theme_set(theme_classic()+
                              theme(axis.ticks.length=unit(.15, "cm"), 
                                    axis.text.x = element_text(angle=45, hjust=1), 
                                    axis.line = element_line(colour = 'black', size = 0.2), 
                                    axis.ticks = element_line(color="black", size=0.2), strip.background =element_rect(fill="black"),
                                    strip.text = element_text(colour = 'white')))
library(ggExtra)
group.colors <- c(CC = "#D55E00", CS = "#E69F00", CM ="#F0E442", CComm = "#56B4E9", CConv = "#0072B2", HV = "#009E73", Flu = "#999999", Sepsis = "#CC79A7", LCC = "#000000")

##load bulk
db <- data.table::fread("/Users/bosun/Documents/GitHub/COMBAT_rep/code/data/All_raw_values_COMBAT_BCR1_VDJ_FINAL.txt")
long <- db[,c(2,3,16:210)] %>%
  reshape2::melt(id.vars=c("RNASeq_sample_ID", "Source"),  na.rm = FALSE) %>%
  filter(Source != "Batch control") %>%
  dplyr::rename(bulk_v = variable, bulk_prop=value)
long$Source <- factor(long$Source, levels=c("HV", "COVID_MILD", "COVID_SEV", "COVID_CRIT", "COVID_HCW_MILD", "Sepsis"))


##load scRNA
Bcell_df <- data.table::fread("/Users/bosun/Documents/GitHub/COMBAT_rep/data/Bcell_df_v001.csv")
##Further clean of data and set factors----
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
Heavy <- Bcell_df[c(Bcell_df$COMBAT_ID_Time %in% samplenames$COMBAT_ID_Time), ]
Heavy <- Heavy %>% 
  filter(!grepl("Flu",Source_abrev)) 
new_source <- Heavy %>% select(Source, Source_abrev) %>% unique
new_source$Source_abrev <- factor(new_source$Source_abrev, levels = c("HV", "CM","CS", "CC","CComm","Flu", "Sepsis", "LCC"))

##filter for those in single cell
long <- db[,c(1,2,3,16:210)] %>%
  reshape2::melt(id.vars=c("Sequencing ID","RNASeq_sample_ID", "Source"),  na.rm = FALSE) %>%
  filter(Source != "Batch control") %>%
  dplyr::rename(bulk_v = variable, bulk_prop=value)
long$Source <- factor(long$Source, levels=c("HV", "COVID_MILD", "COVID_SEV", "COVID_CRIT", "COVID_HCW_MILD", "Sepsis"))

longfilt<- long[c(long$RNASeq_sample_ID %in% Heavy$RNASeq_sample_ID),]
numbers <- longfilt %>% select(Source, RNASeq_sample_ID) %>% unique %>% group_by(Source) %>% add_tally %>% select(n) %>% unique
print(numbers)


subsetlong <- longfilt %>%
  select(RNASeq_sample_ID, Source, `Sequencing ID`)


## Class switched population
#separate into class switched, ighd/ighm unmut, ighd/ighm mut
csw <- longfilt %>%
  filter(grepl("Class_switched", bulk_v)) %>%
  mutate(bulk_v = gsub(".*\\..", "I", bulk_v)) 
sc_csw <- Heavy %>% 
  filter(!grepl("IGHD|IGHM|NA", c_gene)) %>%
  group_by(RNASeq_sample_ID)%>%
  add_tally(name = "sum") %>%
  group_by(RNASeq_sample_ID,v_gene_HC) %>% 
  add_tally(name="sum_p") %>%
  select(RNASeq_sample_ID, Source ,sum,sum_p) %>% 
  unique %>%
  mutate(prop=sum_p/sum) %>%
  arrange(Source,RNASeq_sample_ID) %>%
  select(-sum, -sum_p) %>%
  spread(v_gene_HC, value=prop) 

sc_csw[is.na(sc_csw)] <- 0
sc_csw <- sc_csw %>% 
  reshape2::melt(id.vars = c("RNASeq_sample_ID", "Source"))%>%
  dplyr::rename(sc_csw=variable, sc_prop=value)

combined_csw <- csw %>%
  left_join(sc_csw, by=c("RNASeq_sample_ID","Source"="Source", "bulk_v"="sc_csw"))

combined_csw[is.na(combined_csw)] <- 0

ighv_csw <- combined_csw %>%
  left_join(new_source, by=c("Source")) %>%
  left_join(numbers, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(bulk_prop+0.01, sc_prop+0.01))+
  ggpointdensity::geom_pointdensity(size=1)+
  viridis::scale_color_viridis()+
  theme(legend.position = "none")+
  scale_y_log10()+
  scale_x_log10()+
  xlab("log prop of repertoire by bulk")+
  ylab("log prop of repertoire by single cell")+
  ggtitle("Class switched")
ighv_csw_marginal <- ggMarginal(ighv_csw, type="histogram")
print(ighv_csw_marginal)

boxplot_missed<- combined_csw %>%
  left_join(numbers, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  filter(bulk_prop > 0) %>%
  filter(sc_prop ==0) %>%
  ggplot(aes(bulk_v, bulk_prop))+
  geom_boxplot()+
  facet_wrap(~label)+
  coord_flip()

print(boxplot_missed)

## IGHD/M mutated
dm_mut <- longfilt %>%
  filter(grepl("mutated", bulk_v))%>%
  mutate(bulk_v = gsub(".*\\..", "I", bulk_v)) 

sc_dm_mut <- Heavy %>% 
  filter(grepl("IGHD|IGHM", c_gene)) %>%
  filter(total_mut >4) %>%
  group_by(RNASeq_sample_ID)%>%
  add_tally(name = "sum") %>%
  group_by(RNASeq_sample_ID,v_gene_HC) %>% 
  add_tally(name="sum_p") %>%
  select(RNASeq_sample_ID, Source ,sum,sum_p) %>% 
  unique %>%
  mutate(prop=sum_p/sum) %>%
  arrange(Source,RNASeq_sample_ID) %>%
  select(-sum, -sum_p) %>%
  spread(v_gene_HC, value=prop) 

sc_dm_mut[is.na(sc_dm_mut)] <- 0
sc_dm_mut <- sc_dm_mut %>% 
  reshape2::melt(id.vars = c("RNASeq_sample_ID", "Source"))%>%
  dplyr::rename(sc_dm_mut=variable, sc_prop=value)

combined_dm_mut <- dm_mut %>%
  left_join(sc_dm_mut, by=c("RNASeq_sample_ID","Source"="Source", "bulk_v"="sc_dm_mut"))

combined_dm_mut[is.na(combined_dm_mut)] <- 0
ighv_dm_mut <- combined_dm_mut %>%
  left_join(new_source, by=c("Source")) %>%
  left_join(numbers, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(bulk_prop+0.01, sc_prop+0.01))+
  ggpointdensity::geom_pointdensity(size=1)+
  viridis::scale_color_viridis()+
  theme(legend.position = "none")+
  scale_y_log10()+
  scale_x_log10()+
  xlab("log prop of repertoire by bulk")+
  ylab("log prop of repertoire by single cell")+
  ggtitle("IGHD/M mutated")
ighv_dm_mut_marginal <- ggMarginal(ighv_dm_mut, type="histogram")

boxplot_missed_dm_mut<- combined_dm_mut %>%
  left_join(numbers, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  filter(bulk_prop > 0) %>%
  filter(sc_prop ==0) %>%
  ggplot(aes(bulk_v, bulk_prop))+
  geom_boxplot()+
  facet_wrap(~label)+
  coord_flip()

print(ighv_dm_mut)

print(boxplot_missed_dm_mut)

## IGHD/M unmutated
dm_unmut <- longfilt %>%
  filter(grepl("unmutated", bulk_v))%>%
  mutate(bulk_v = gsub(".*\\..", "I", bulk_v)) 

sc_dm_unmut <- Heavy %>% 
  filter(grepl("IGHD|IGHM", c_gene)) %>%
  filter(total_mut <=4) %>%
  group_by(RNASeq_sample_ID)%>%
  add_tally(name = "sum") %>%
  group_by(RNASeq_sample_ID,v_gene_HC) %>% 
  add_tally(name="sum_p") %>%
  select(RNASeq_sample_ID, Source ,sum,sum_p) %>% 
  unique %>%
  mutate(prop=sum_p/sum) %>%
  arrange(Source,RNASeq_sample_ID) %>%
  select(-sum, -sum_p) %>%
  spread(v_gene_HC, value=prop) 

sc_dm_unmut[is.na(sc_dm_unmut)] <- 0
sc_dm_unmut <- sc_dm_unmut %>% 
  reshape2::melt(id.vars = c("RNASeq_sample_ID", "Source"))%>%
  dplyr::rename(sc_dm_unmut=variable, sc_prop=value)

combined_dm_unmut <- dm_unmut %>%
  left_join(sc_dm_unmut, by=c("RNASeq_sample_ID","Source"="Source", "bulk_v"="sc_dm_unmut"))

combined_dm_unmut[is.na(combined_dm_unmut)] <- 0
ighv_dm_unmut <- combined_dm_unmut %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(bulk_prop+0.01, sc_prop+0.01))+
  ggpointdensity::geom_pointdensity(size=1)+
  viridis::scale_color_viridis()+
  theme(legend.position = "none")+  
  scale_y_log10()+
  scale_x_log10()+
  xlab("log prop of repertoire by bulk")+
  ylab("log prop of repertoire by single cell")+
  ggtitle("IGHD/M unmutated")
ighv_dm_unmut_marginal <- ggMarginal(ighv_dm_unmut, type="histogram")

boxplot_missed_dm_unmut<- combined_dm_unmut %>%
  left_join(numbers, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  filter(bulk_prop > 0) %>%
  filter(sc_prop ==0) %>%
  ggplot(aes(bulk_v, bulk_prop))+
  geom_boxplot()+
  facet_wrap(~label)+
  coord_flip()

print(ighv_dm_unmut_marginal)

print(boxplot_missed_dm_unmut)

#Compare single cell vs bulk
sc <- combined_dm_unmut %>%
  filter(bulk_prop==0) %>%
    filter(sc_prop!=0) %>%
   left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(sc_prop*100, fill=Source_abrev))+
  geom_histogram()+
  facet_wrap(~Source_abrev, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_x_log10()+ 
  ggtitle("IGHD/M unmutated - missed by bulk")+
  theme(legend.position = "none")


bulk <- combined_dm_unmut %>%
  filter(sc_prop==0) %>%
      filter(bulk_prop!=0) %>%
   left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(bulk_prop, fill=Source_abrev))+
  geom_histogram()+
  facet_wrap(~Source_abrev, nrow=1)+
scale_fill_manual(values=group.colors)+
  scale_x_log10()+
  scale_x_log10()+ ggtitle("IGHD/M unmutated - missed by single cell")+
  theme(legend.position = "none")


p_unmut<- patchwork::wrap_plots(sc/bulk) 


bulk_missed_dm_unmut <- combined_dm_unmut %>%
  filter(bulk_prop==0) %>%
  filter(sc_prop!=0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  select(RNASeq_sample_ID,Source_abrev, bulk_v) %>%
  unique() %>%
  group_by(RNASeq_sample_ID) %>%
  add_tally() %>%
  select(Source_abrev, n) %>%
  unique %>%
  ggplot(aes(Source_abrev, n, fill=Source_abrev))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=1)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  ylab("Number of V genes missed")+
  ggtitle("Missed by bulk")+
  ylim(0,45)+
  theme(legend.position = "none")
sc_missed_dm_unmut <- combined_dm_unmut %>%
  filter(bulk_prop!=0) %>%
  filter(sc_prop==0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  select(RNASeq_sample_ID,Source_abrev, bulk_v) %>%
  unique() %>%
  group_by(RNASeq_sample_ID) %>%
  add_tally() %>%
  select(Source_abrev, n) %>%
  unique %>%
  ggplot(aes(Source_abrev, n, fill=Source_abrev))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=1)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  ylab("Number of V genes missed")+
  ggtitle("Missed by single cell")+
  ylim(0,45)+
  theme(legend.position = "none")

dm_unmut_comb<- bulk_missed_dm_unmut|sc_missed_dm_unmut

sc2 <- combined_dm_mut %>%
  filter(bulk_prop==0) %>%
  filter(sc_prop!=0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(sc_prop*100, fill=Source_abrev))+
  geom_histogram()+
  facet_wrap(~Source_abrev, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_x_log10()+ ggtitle("IGHD/M mutated - missed by bulk")+
  theme(legend.position = "none")

bulk2 <- combined_dm_mut %>%
  filter(sc_prop==0) %>%
  filter(bulk_prop!=0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(bulk_prop, fill=Source_abrev))+
  geom_histogram()+
  facet_wrap(~Source_abrev, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_x_log10()+
  scale_x_log10()+ ggtitle("IGHD/M mutated - missed by single cell")+
  theme(legend.position = "none")

p_mut<- patchwork::wrap_plots(sc2/bulk2) 

bulk_missed_dm_mut <- combined_dm_mut %>%
  filter(bulk_prop==0) %>%
  filter(sc_prop!=0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  select(RNASeq_sample_ID,Source_abrev, bulk_v) %>%
  unique() %>%
  group_by(RNASeq_sample_ID) %>%
  add_tally() %>%
  select(Source_abrev, n) %>%
  unique %>%
  ggplot(aes(Source_abrev, n, fill=Source_abrev))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=1)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  ylab("Number of V genes missed")+
  ggtitle("Missed by bulk")+
  ylim(0,45)+
  theme(legend.position = "none")
sc_missed_dm_mut <- combined_dm_mut %>%
  filter(bulk_prop!=0) %>%
  filter(sc_prop==0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  select(RNASeq_sample_ID,Source_abrev, bulk_v) %>%
  unique() %>%
  group_by(RNASeq_sample_ID) %>%
  add_tally() %>%
  select(Source_abrev, n) %>%
  unique %>%
  ggplot(aes(Source_abrev, n, fill=Source_abrev))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=1)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  ylab("Number of V genes missed")+
  ggtitle("Missed by single cell")+
  ylim(0,45)+
  theme(legend.position = "none")

dm_mut_comb<- bulk_missed_dm_mut|sc_missed_dm_mut

sc3 <- combined_csw %>%
  filter(bulk_prop==0) %>%
  filter(sc_prop!=0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(sc_prop*100, fill=Source_abrev))+
  geom_histogram()+
  facet_wrap(~Source_abrev, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_x_log10()+ ggtitle("Class switched - missed by bulk")+
  theme(legend.position = "none")


bulk3 <- combined_csw %>%
  filter(sc_prop==0) %>%
  filter(bulk_prop!=0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  mutate(label = paste(Source, n, sep =" n=")) %>%
  ggplot(aes(bulk_prop, fill=Source_abrev))+
  geom_histogram()+
  facet_wrap(~Source_abrev, nrow=1)+
  scale_fill_manual(values=group.colors)+
  scale_x_log10()+ ggtitle("Class switched - missed by single cell")+
  theme(legend.position = "none")

p_csw<- patchwork::wrap_plots(sc3/bulk3) 

bulk_missed_csw <- combined_csw %>%
  filter(bulk_prop==0) %>%
  filter(sc_prop!=0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  select(RNASeq_sample_ID,Source_abrev, bulk_v) %>%
  unique() %>%
  group_by(RNASeq_sample_ID) %>%
  add_tally() %>%
  select(Source_abrev, n) %>%
  unique %>%
  ggplot(aes(Source_abrev, n, fill=Source_abrev))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=1)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  ylab("Number of V genes missed")+
  ggtitle("Missed by Bulk")+
  ylim(0,45)+
  theme(legend.position = "none")
sc_missed_csw <- combined_csw %>%
  filter(bulk_prop!=0) %>%
  filter(sc_prop==0) %>%
  left_join(numbers, by=c("Source")) %>%
  left_join(new_source, by=c("Source")) %>%
  select(RNASeq_sample_ID,Source_abrev, bulk_v) %>%
  unique() %>%
  group_by(RNASeq_sample_ID) %>%
  add_tally() %>%
  select(Source_abrev, n) %>%
  unique %>%
  ggplot(aes(Source_abrev, n, fill=Source_abrev))+
  geom_boxplot(aes(fill=Source_abrev),width=0.5, alpha=0.3, outlier.shape = NA, coef=0, color=NA) +
  geom_jitter(aes(color=Source_abrev),width=0.1, size=1)+
  stat_summary(aes(color=Source_abrev), fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.6, alpha=90)+
  scale_fill_manual(values=group.colors)+
  scale_color_manual(values=group.colors)+
  ylab("Number of V genes missed")+
  ggtitle("Missed by single cell")+
  ylim(0,45)+
  theme(legend.position = "none")

csw_comb<- bulk_missed_csw|sc_missed_csw

combined_sc_bulk <- patchwork::wrap_plots(p_unmut, p_mut, p_csw)

combined_sc_bulk2 <- patchwork::wrap_plots(dm_unmut_comb, dm_mut_comb, csw_comb)

combined_sc_bulk|combined_sc_bulk2

overall_plot <- patchwork::wrap_plots(
  patchwork::wrap_plots(ighv_dm_unmut_marginal,dm_unmut_comb, widths = c(0.8, 1)),
    patchwork::wrap_plots(ighv_dm_mut_marginal,dm_mut_comb, widths=c(0.8,1)),
  patchwork::wrap_plots(ighv_csw_marginal,csw_comb, widths=c(0.8,1)), nrow=3
)
overall_plot

pdf("COMBAT_figure_plots/ighv_qc_plot.pdf", heigh=10, width=8)
overall_plot
dev.off()


##stats----
new_labels <- Heavy %>% select(Source, Source_abrev) %>% unique
sc_dm_unmut <- sc_dm_unmut %>% 
  left_join(new_labels, by=c("Source"))
sc_dm_unmut_dunn <- sc_dm_unmut %>%
  group_by(sc_dm_unmut) %>%
  rstatix::dunn_test(sc_prop ~ Source_abrev, p.adjust.method = "fdr")

sc_dm_unmut_list <- sc_dm_unmut %>%
  group_by(sc_dm_unmut) %>%
  group_split() %>%
  setNames(unique(sc_dm_unmut$sc_dm_unmut))

sc_dm_unmut_long_mean_logfc <- lapply(sc_dm_unmut_list, function(x){
  pairwise.fc(x$sc_prop, x$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
    data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
    mutate(sc_dm_unmut=unique(x$sc_dm_unmut)) %>%
    mutate(condition=paste(condition1, sc_dm_unmut, sep="-")) %>%
    separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
})

sc_dm_unmut_long_mean_logfc <- data.table::rbindlist(sc_dm_unmut_long_mean_logfc)

sc_dm_unmut_dunn <- sc_dm_unmut_dunn %>% 
  left_join(sc_dm_unmut_long_mean_logfc, by=c("group1", "group2", "sc_dm_unmut")) 

sc_dm_unmut_dunn_plot <- sc_dm_unmut_dunn %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.3)+
  geom_point(shape=16, alpha=0.6)+
  ggrepel::geom_text_repel(data=filter(sc_dm_unmut_dunn, p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  theme(legend.position = "none")+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))


sc_dm_unmut_dunn %>%
  filter(sc_dm_unmut == "IGHV3-33") %>%
  select(group1, group2, statistic, fold_change, p.adj.signif)

pdf("COMBAT_figure_plots/ighv_igdm_unmutated_volcano.pdf", height=7.6, width=6)
sc_dm_unmut_dunn_plot
dev.off()

#sanity check plot
sc_dm_unmut %>%
  filter(sc_dm_unmut == "IGHV3-33") %>%
  ggplot(aes(Source_abrev, sc_prop))+
  geom_boxplot()

##stats for IGHD/M mutated -bulk

#single cell
sc_dm_mut %>%
  group_by(sc_dm_mut) %>%
  rstatix::dunn_test(sc_prop~Source, p.adjust.method = "fdr") %>%
  filter(p.adj.signif != "ns")

#bulk
dm_mut <- dm_mut %>% 
  left_join(new_labels, by=c("Source"))
dm_mutlist <-  split(dm_mut, dm_mut$bulk_v)

dm_mut_dunn <- lapply(dm_mutlist, function(x){
  if(sum(x$bulk_prop) != 0){
    x %>%
      rstatix::dunn_test(bulk_prop~Source_abrev, p.adjust.method = "fdr") %>%
      mutate(bulk_v=unique(x$bulk_v))
  }else{
    NULL
  }
})
 
dm_mut_dunn[sapply(dm_mut_dunn, is.null)] <- NULL

dm_mut_dunn <- data.table::rbindlist(dm_mut_dunn)

dm_mut_list <- dm_mut %>%
  group_by(bulk_v) %>%
  group_split() %>%
  setNames(unique(dm_mut$bulk_v))

dm_mut_long_mean_logfc <- lapply(dm_mut_list, function(x){
  pairwise.fc(x$bulk_prop, x$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
    data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
    mutate(bulk_v=unique(x$bulk_v)) %>%
    mutate(condition=paste(condition1, bulk_v, sep="-")) %>%
    separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
})

dm_mut_long_mean_logfc <- data.table::rbindlist(dm_mut_long_mean_logfc)

dm_mut_dunn <- dm_mut_dunn %>% 
  left_join(dm_mut_long_mean_logfc, by=c("group1", "group2", "bulk_v")) 

dm_mut_dunn_plot <- dm_mut_dunn %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.3)+
  geom_point(shape=16, alpha=0.6)+
  ggrepel::geom_text_repel(data=filter(dm_mut_dunn, p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  theme(legend.position = "none")+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))

pdf("COMBAT_figure_plots/ighv_igdm_mutated_volcano.pdf", height=7.6, width=6)
dm_mut_dunn_plot
dev.off()

##CSW stats ----
csw <- csw %>% 
  left_join(new_labels, by=c("Source"))

csw_mutlist <-  split(csw, csw$bulk_v)

csw_dunn <- lapply(csw_mutlist, function(x){
  if(sum(x$bulk_prop) != 0){
    x %>%
      rstatix::dunn_test(bulk_prop~Source_abrev, p.adjust.method = "fdr") %>%
      mutate(bulk_v=unique(x$bulk_v))
  }else{
    NULL
  }
})

csw_dunn[sapply(csw_dunn, is.null)] <- NULL

csw_dunn <- data.table::rbindlist(csw_dunn)

csw_list <- csw %>%
  group_by(bulk_v) %>%
  group_split() %>%
  setNames(unique(csw$bulk_v))
library(MKmisc)

csw_long_mean_logfc <- lapply(csw_list, function(x){
  pairwise.fc(x$bulk_prop, x$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
    data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
    mutate(bulk_v=unique(x$bulk_v)) %>%
    mutate(condition=paste(condition1, bulk_v, sep="-")) %>%
    separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
})

csw_long_mean_logfc <- data.table::rbindlist(csw_long_mean_logfc)

csw_dunn <- csw_dunn %>% 
  left_join(csw_long_mean_logfc, by=c("group1", "group2", "bulk_v")) 

csw_dunn_plot <- csw_dunn %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.3)+
  geom_point(shape=16, alpha=0.6)+
  ggrepel::geom_text_repel(data=filter(csw_dunn, p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  theme(legend.position = "none")+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))

pdf("COMBAT_figure_plots/ighv_igcswated_volcano.pdf", height=7.6, width=6)
csw_dunn_plot
dev.off()

##Attempt to combine plots
sc_dm_umut_compat <- sc_dm_unmut_dunn %>%
  dplyr::rename(bulk_v=sc_dm_unmut)
combined <- rbind(sc_dm_umut_compat %>% mutate(label="IGHD/M unmutated (single cell)"), 
      dm_mut_dunn%>% mutate(label="IGHD/M mutated (RNA)"), 
      csw_dunn%>% mutate(label="Class switched (RNA)"))

#filter out sepsis
combined$label <- factor(combined$label, levels=c("IGHD/M unmutated (single cell)", "IGHD/M mutated (RNA)", "Class switched (RNA)"))

covidvsHV <- combined %>%
  filter(!grepl("Sepsis", group1)) %>%
  filter(!grepl("Sepsis", group2)) %>%
  filter(grepl("HV", group2)) %>%
  mutate(overall="COVID-19 vs. HV") 

condcomb <- rbind(covidvssepsis, covidvsHV)
condcomb$label <- factor(condcomb$label, levels=c("IGHD/M unmutated (single cell)", "IGHD/M mutated (RNA)", "Class switched (RNA)"))

covidvsHV_plot <- covidvsHV %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.3)+
  geom_point(shape=16, alpha=0.6)+
  ggrepel::geom_text_repel(data=filter(covidvsHV , p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  theme(legend.position = "none")+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))+
  facet_wrap(~label, scales="free_y")+
  ggtitle(paste(unique(covidvsHV$overall)))


pdf("COMBAT_figure_plots/all_ighv_covidvsHV_volcano.pdf", height=5.5, width=13)
covidvsHV_plot
dev.off()

covidvssepsis <- combined %>%
  filter(!grepl("HV", group1)) %>%
  filter(!grepl("HV", group2)) %>%
  filter(grepl("Sepsis", group2)) %>%
  mutate(overall="COVID-19 vs. Sepsis")

covidvssepsis_plot <- covidvssepsis %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.3)+
  geom_point(shape=16, alpha=0.6)+
  ggrepel::geom_text_repel(data=filter(covidvssepsis , p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  theme(legend.position = "none")+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))+
  facet_wrap(~label, scales="free_y")+
  ggtitle(paste(unique(covidvssepsis$overall)))

pdf("COMBAT_figure_plots/all_ighv_covidvsSepsis_volcano.pdf", height=5.5, width=13)
covidvssepsis_plot
dev.off()

covidvscovid <- combined %>%
  filter(!grepl("HV", group1)) %>%
  filter(!grepl("HV", group2)) %>%
  filter(!grepl("Sepsis", group1)) %>%
  filter(!grepl("Sepsis", group2)) %>%
  mutate(overall="CComm vs CM vs CS vs CC")

covidvscovid_plot <- covidvscovid %>%
  mutate(col= ifelse(p.adj<0.05 & fold_change >1, "red", ifelse(p.adj<0.05 & fold_change <1, "blue", "gray"))) %>%
  ggplot(aes(log2(fold_change), -log10(p.adj), color=col))+
  geom_point(shape=19, color="black", alpha=0.3)+
  geom_point(shape=16, alpha=0.6)+
  ggrepel::geom_text_repel(data=filter(covidvscovid , p.adj<0.05), 
                           aes(log2(fold_change), -log10(p.adj), label=condition), color="black", size=3, max.overlaps = Inf, 
                           segment.size	=0.3)+
  geom_vline(xintercept = 0, linetype="dashed", alpha=0.2)+
  scale_color_manual(values=c("#2171b5", "gray", "#cb181d"))+
  xlab("avg_logFC of repertoire proportion")+
  theme(legend.position = "none")+
  theme(legend.position = "none", axis.text.x=element_text(angle=0, hjust=0.5))+
  facet_wrap(~label, scales="free_y")+
  ggtitle(paste(unique(covidvscovid$overall)))

pdf("COMBAT_figure_plots/all_ighv_covidvscovid_volcano.pdf", height=5.5, width=13)
covidvscovid_plot
dev.off()

covidvscovid %>%
  filter(p.adj.signif != "ns")

##Bulk VJ-----
vj_bulk <- as.data.frame(data.table::fread("data/All_Mutations_VJ_genes.txt"))
vj_bulk <- vj_bulk %>%
  filter(!grepl("TCR", `#sample`))
new_bulk<- subsetlong %>% unique %>%
  left_join(vj_bulk, by=c("Sequencing ID"="#sample"))

new_bulk <- new_bulk %>% 
  left_join(new_labels, by=c("Source"))

bulkisoVJ <- new_bulk %>% 
  filter(grepl("IGHD/M_|Class_switched", isotype)) 

bulkisoVJ <- bulkisoVJ %>%
  mutate(vj_pair= paste(v,j, sep="|"))

bulkisoVJ$n_unique_BCRs <- as.numeric(bulkisoVJ$n_unique_BCRs)


total_BCR <- bulkisoVJ %>%
  group_by(RNASeq_sample_ID, isotype) %>%
  dplyr::summarise(n_total_bcrs = sum(n_unique_BCRs))

bulkisoVJ <- bulkisoVJ %>%
  left_join(total_BCR, by=c("RNASeq_sample_ID", "isotype"))

bulkisoVJ <- bulkisoVJ %>%
  mutate(prop_of_total= n_unique_BCRs/n_total_bcrs)

unique(bulkisoVJ$RNASeq_sample_ID) %in% Heavy$RNASeq_sample_ID %>% summary
unique(longfilt$RNASeq_sample_ID) %in% Heavy$RNASeq_sample_ID %>% summary
(Heavy$RNASeq_sample_ID %>%unique) %in%longfilt$RNASeq_sample_ID %>% summary

bulkisoVJ <- bulkisoVJ[c(bulkisoVJ$RNASeq_sample_ID %in% Heavy$RNASeq_sample_ID),]
bulkisovjlist <- split(bulkisoVJ, bulkisoVJ$isotype)


##stats

dunnbulkvj <- lapply(bulkisovjlist, function(x){
  y <-  split(x, x$vj_pair)
  
  dunn <- lapply(y, function(z){
    if(sum(z$prop_of_total) != 0 & length(unique(z$Source_abrev))>1){
      z %>%
        ungroup %>% 
        rstatix::dunn_test(prop_of_total~Source_abrev, p.adjust.method = "fdr") %>%
        mutate(vj_pair=unique(z$vj_pair))
    }else{
      NULL
    }
  })
  
  dunn[sapply(dunn, is.null)] <- NULL
  dunn <- data.table::rbindlist(dunn)
 
   long_mean_logfc <- lapply(y, function(v){
    if(sum(v$prop_of_total) != 0 & length(unique(v$Source_abrev))>1){
      pairwise.fc(v$prop_of_total, v$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
        data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
        mutate(vj_pair=unique(v$vj_pair)) %>%
        mutate(condition=paste(condition1, vj_pair, sep="-")) %>%
        separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
    }else{
      NULL
    }
  })
  
  long_mean_logfc <- data.table::rbindlist(long_mean_logfc)
  
  dunn <- dunn %>% 
    left_join(long_mean_logfc, by=c("group1", "group2", "vj_pair")) 
  
  return(dunn)
})

bulkvjnames<- names(dunnbulkvj)

vjbulk_plotlist <- lapply(bulkvjnames, function(x){
  plot <- dunnbulkvj[[x]]  %>%
      mutate(group1vs=paste(group1, "vs.", sep=" ")) %>%
      filter(!grepl("ns", p.adj.signif)) %>%
      filter(!grepl("Inf", fold_change)) %>%
      filter(fold_change != 0) %>%
      ggplot(aes(vj_pair, group2))+
      geom_point(aes(size=log2(fold_change),color=log2(fold_change)), alpha=0.8, shape=19)+
      geom_text(aes(label=p.adj.signif), vjust=0.8)+
      facet_wrap(~group1vs, nrow=1)+
      ggsci::scale_color_gsea()+
      coord_flip()+
      ylab("")+
      xlab("IGHV-J pairing")+
    ggtitle(paste0(x))+
    theme(legend.position = "bottom")
  
  return(plot)
})

##single cell
sc_dm_unmut_vj <- Heavy %>% 
  filter(grepl("IGHD|IGHM", c_gene)) %>%
  filter(total_mut <=4) %>%
  group_by(RNASeq_sample_ID)%>%
  add_tally(name = "sum") %>%
  group_by(RNASeq_sample_ID,VJ_HC) %>% 
  add_tally(name="sum_p") %>%
  select(RNASeq_sample_ID, Source ,sum,sum_p) %>% 
  unique %>%
  mutate(prop=sum_p/sum) %>%
  arrange(Source,RNASeq_sample_ID) %>%
  select(-sum, -sum_p) %>%
  spread(VJ_HC, value=prop) 

sc_dm_unmut_vj[is.na(sc_dm_unmut_vj)] <- 0
sc_dm_unmut_vj <- sc_dm_unmut_vj %>% 
  reshape2::melt(id.vars = c("RNASeq_sample_ID", "Source"))%>%
  dplyr::rename(vj_pair=variable, sc_prop=value)

sc_dm_unmut_vj <- sc_dm_unmut_vj %>%
  left_join(new_labels, by=c("Source"))

dunn_sc_unmutvj <- sc_dm_unmut_vj %>%
  group_by(vj_pair) %>%
  rstatix::dunn_test(sc_prop ~ Source_abrev, p.adjust.method = "fdr")

scunmutvj_list <- sc_dm_unmut_vj %>%
  group_by(vj_pair) %>%
  group_split() %>%
  setNames(unique(sc_dm_unmut_vj$vj_pair))

sc_unmutvj_fold <- lapply(scunmutvj_list, function(x){
  pairwise.fc(x$sc_prop, x$Source_abrev, ave = mean, log = FALSE,base=2, mod.fc = FALSE) %>%
    data.frame(condition1=names(.), fold_change=., row.names=NULL) %>%
    mutate(vj_pair=unique(x$vj_pair)) %>%
    mutate(condition=paste(condition1, vj_pair, sep="-")) %>%
    separate(condition1, c("group1", "group2"), sep="[[:space:]]vs[[:space:]]")
})

dunn_sc_unmutvj <- dunn_sc_unmutvj %>% 
  left_join(sc_unmutvj_fold, by=c("group1", "group2", "vj_pair")) 

sc_vj_unmut_plot <- dunn_sc_unmutvj %>%
  mutate(group1vs=paste(group1, "vs.", sep=" ")) %>%
  filter(!grepl("ns", p.adj.signif)) %>%
  filter(!grepl("Inf", fold_change)) %>%
  filter(fold_change != 0) %>%
  ggplot(aes(vj_pair, group2))+
  geom_point(aes(size=log2(fold_change),color=log2(fold_change)), alpha=0.8, shape=19)+
  geom_text(aes(label=p.adj.signif), vjust=0.8)+
  facet_wrap(~group1vs, nrow=1)+
  ggsci::scale_color_gsea()+
  coord_flip()+
  ylab("")+
  xlab("IGHV-J pairing")+
  ggtitle("IGHD/M unmutated")+
  theme(legend.position = "bottom")

vjbulk_plotlist[[3]] <- sc_vj_unmut_plot

pdf("COMBAT_figure_plots/VJ_plots.pdf", height=16, width=9)
patchwork::wrap_plots(vjbulk_plotlist, ncol=1)
dev.off()




