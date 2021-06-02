##Light chain swapping script
l_locus <- data.table::fread("/Users/bosun/Documents/GitHub/COMBAT_rep/code/data/light_chain_locus_imgt.csv")

l_locus <- l_locus %>%
  separate(IMGT_position, c("locus_position_from", "locus_position_to"))
l_locus$locus_position_from <- as.numeric(l_locus$locus_position_from)
l_locus$locus_position_to <- as.numeric(l_locus$locus_position_to)



Light$genomic_distance <- l_locus[match(Light$j_gene_LC, l_locus$Gene)]$locus_position_from-l_locus[match(Light$v_gene_LC, l_locus$Gene)]$locus_position_from

length(Light$v_gene_LC)
length(Light$j_gene_LC)

length(l_locus[match(Light$v_gene_LC, l_locus$Gene)]$locus_position_from)
length(l_locus[match(Light$j_gene_LC, l_locus$Gene)]$locus_position_from)

Light$Source <- factor(Light$Source, levels=c("HV", "COVID_MILD", "COVID_SEV", "COVID_CRIT", "COVID_HCW_MILD", "Sepsis", "Flu"))

Light <- Light[c(Light$barcode_id %in% Heavy$barcode_id),]

Light <- Light %>%
  select(-pseudobulk) %>%
  left_join(Heavy %>% select(barcode_id, pseudobulk), by=c("barcode_id"))

k_plot <- aggregate(genomic_distance ~ COMBAT_ID_Time+Source+pseudobulk, data=Light %>% filter(grepl("IGK", c_gene)), mean) %>%
  ggplot(aes(Source, genomic_distance))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~pseudobulk, scales="free_y")+
  scale_y_log10()+
  ggtitle("Kappa genomic distance")

k_tukey <- aggregate(genomic_distance ~ COMBAT_ID_Time+Source+pseudobulk, data=Light %>% filter(grepl("IGK", c_gene)), mean) %>%
  group_by(pseudobulk) %>%
  rstatix::tukey_hsd(genomic_distance ~ Source)
k_tukey %>% filter(p.adj.signif != "ns")
l_plot <- aggregate(genomic_distance ~ COMBAT_ID_Time+Source+pseudobulk, data=Light %>% filter(grepl("IGL", c_gene)), mean) %>%
  ggplot(aes(Source, genomic_distance))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~pseudobulk, scales="free_y")+
  ggtitle("Lambda genomic distance")

l_tukey <- aggregate(genomic_distance ~ COMBAT_ID_Time+Source+pseudobulk, data=Light %>% filter(grepl("IGL", c_gene)), mean) %>%
  group_by(pseudobulk) %>%
  rstatix::tukey_hsd(genomic_distance ~ Source) %>%
  add_xy_position(x="Source")

l_tukey %>% filter(p.adj.signif != "ns")


l_plot <- l_plot+ ggpubr::stat_pvalue_manual(l_tukey, hide.ns = TRUE  , label = "p.adj.signif", tip.length = 0.01 , step.increase = 0.1,)


patchwork::wrap_plots(k_plot, l_plot)

Light %>%
  ggplot(aes(genomic_distance))+
  geom_histogram()+
  facet_wrap(~COMBAT_ID_Time)


