##Script to generate meta data from BCR repertoire to guide clustering annotations
##Bo Sun
##Bashford-Rogers Group
##Date: 2020-08-27

require(tidyverse)
require(pathwork)
require(ggsci)
require(RColorBrewer)

#Load clonally assigned heavy chains----
#Clones were assigned using BcellR- unpublished single cell repertoire analysis tool (Sun et al)
#In brief, VJ nucleotide sequence was fuzzy clustered using cdhit, then normalised levenhsteins distances were calculated pairwise across all CDR3s
#Finally, for speed, a density model based optimization was applied to find the local minimum of the binomial distribution of CDR3 distances. 
clones <- readRDS("/well/combat/users/vkh192/repertoire/data/new_clones.rds")

#Naive B cells are defined as those with no mutations, although V gene rearrangement can mean alignment tools may assign 1-2 mutations at the junction region
#Here we allow for up to 2 mutations from germline
naive <- clones$Heavy %>%
  filter(grepl("IGHM|IGHD", c_gene)) %>%
  filter(total_mut < 2) %>% 
  select(barcode_seq, c_gene, COMBAT_ID, COMBAT_ID_Time, name, subset)%>%
  dplyr::rename(barcode_id = barcode_seq) %>%
  mutate(new_name = "naive.B")

#Activated B cells were defined as near germline B cells that had either class switched or obtained more than n=2 mutations. 
#These can also be considered to be extrafollicular B cells of DN2 phenotype if they have the following surface phenotype: 
#(IgD-CD27-CD38-CD24-CD21-CXCR5-Tbet+CD11c+FcRL5+SLAMF7+)
activatedB <- clones$Heavy %>%
  filter(grepl("naive", name))%>%
  filter(total_mut > 2) %>% 
  select(barcode_seq, c_gene, COMBAT_ID, COMBAT_ID_Time, name, subset) %>%
  dplyr::rename(barcode_id = barcode_seq) %>%
  mutate(new_name = "activated.B")

#Plasma cells are quite readily defined within the GEX dataset, so these annotations are really to support the immunoglobulin isotype classification
PC <- clones$Heavy %>% filter(major_cell_type=="PC") %>% 
  select(barcode_seq, c_gene, COMBAT_ID, COMBAT_ID_Time, name, subset) %>%
  dplyr::rename(barcode_id = barcode_seq) %>%
  mutate(iso=gsub("\\*.*", "", c_gene)) %>% #Here I remove alleles, mostly to reduce complexity and also due to low coverage of the reverse primers to make accurate calls
  mutate(new_name = paste(iso, "PC", sep=".")) %>% select(-iso) #new annotation based on BCR repertoire consisting of isotype.PC e.g. IGHA1.PC
  
#All barcodes with c_gene and mutation annotations
BCRmeta <- clones$Heavy %>% 
  select(barcode_seq, total_mut, c_gene, COMBAT_ID, COMBAT_ID_Time, name, subset) %>%
  dplyr::rename(barcode_id = barcode_seq) 


#Write csv
write.csv(naive, "/well/combat/users/vkh192/repertoire/data/VDJ_naive.csv")
write.csv(activatedB, "/well/combat/users/vkh192/repertoire/data/VDJ_activated.csv")
write.csv(PC, "/well/combat/users/vkh192/repertoire/data/VDJ_PC_isotype.csv")
write.csv(BCRmeta, "/well/combat/users/vkh192/repertoire/data/BCRmeta.csv")

##Plots to aid annotation ----
#Load umap
umap <- data.table::fread("/well/combat/datasets/CBD-10X-00004/cell_subsets/d_bcell_plasma/umap.tsv.gz")

#Merge umap
Heavy2 <- Heavy %>% left_join(umap, by=c("barcode_seq"="barcode"))

#Plots for mutations overlaid on umap
p1<- Heavy %>%
  ggplot(aes(UMAP_1, UMAP_2, color=total_mut))+
  geom_point(alpha=0.7, size=0.7)+
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))
p2<- Heavy %>%
  ggplot(aes(UMAP_1, UMAP_2, color=total_mut<2))+
  geom_point(alpha=0.7, size=0.7)+
  scale_color_manual(values=c("gray", "red"))

png("/well/combat/users/vkh192/repertoire/out/umap_mut.png",type=c("cairo"), height = 6, width = 12, units = 'in', res = 400)
p1|p2
dev.off()

#Plots of B cell isotype usage
png("/well/combat/users/vkh192/repertoire/out/umap_iso.png",type=c("cairo"), height = 6, width = 8, units = 'in', res = 400)
Heavy %>%
  mutate(c_gene=gsub("\\*.*", "", c_gene)) %>%
  ggplot(aes(UMAP_1, UMAP_2, fill=c_gene))+
  geom_point(alpha=0.6, size=0.4, shape=21, color="darkgray")+
  ggsci::scale_fill_d3()+
  facet_wrap(~c_gene)
dev.off()

#Plot of c_gene isotype use colored by umis to visualise plasma cells/activated cells with high umis or
#to visualise resting cells which should have low umis
library(ggpointdensity)
#resample so density measure is evenly distributed per bootstrap
Combined <- bind_rows(replicate(100, Heavy2 %>% group_by(COMBAT_ID_Time) %>% sample_n(30, replace=TRUE), simplify=F), .id="COMBAT_ID_Time")

png("/well/combat/users/vkh192/repertoire/out/umap_iso_umis.png",type=c("cairo"), height = 10, width = 11, units = 'in', res = 400)
Combined %>%
  filter(!grepl("None|IGHE", c_gene)) %>%
  mutate(c_gene=gsub("\\*.*", "", c_gene)) %>%
  ggplot(aes(UMAP_1, UMAP_2))+
  geom_pointdensity(alpha=0.6, size=0.4)+
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  facet_grid(vars(c_gene), vars(Source))
dev.off()

png("/well/combat/users/vkh192/repertoire/out/umap_mut2.png",type=c("cairo"), height = 10, width = 11, units = 'in', res = 400)
Heavy2 %>%
  filter(!grepl("None|IGHE", c_gene)) %>%
  mutate(c_gene=gsub("\\*.*", "", c_gene)) %>%
  ggplot(aes(UMAP_1, UMAP_2))+
  geom_point(shape=19, color="black", size=0.5)+
  geom_point(shape=16, aes(color=total_mut), alpha=0.9, size=0.5)+
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  facet_grid(vars(c_gene), vars(Source))
dev.off()

png("/well/combat/users/vkh192/repertoire/out/umap_cdrl.png",type=c("cairo"), height = 10, width = 11, units = 'in', res = 400)
Heavy2 %>%
  filter(!grepl("None|IGHE", c_gene)) %>%
  mutate(c_gene=gsub("\\*.*", "", c_gene)) %>%
  ggplot(aes(UMAP_1, UMAP_2))+
  geom_point(shape=19, color="black", size=0.5)+
  geom_point(shape=16, aes(color=junction_length_HC), alpha=0.9, size=0.5)+
  scale_color_distiller(palette = 'RdYlBu')+
  facet_wrap(~Source)
dev.off()

png("/well/combat/users/vkh192/repertoire/out/umap_clonality.png",type=c("cairo"), height = 10, width = 11, units = 'in', res = 400)
Heavy2 %>%
  ggplot(aes(UMAP_1, UMAP_2))+
  geom_point(shape=19, color="black", size=0.1)+
  geom_point(shape=16, color="gray", alpha=0.5, size=0.5)+
  stat_density2d(data=Heavy2 %>% filter(clonal_abundance>5),aes(alpha=..level.., color=..level..), size=3, bins=10, contour_var = "ndensity")+
  scale_color_distiller(palette = 'RdYlBu')+
  facet_wrap(~Source)
dev.off()

png("/well/combat/users/vkh192/repertoire/out/umap_clonality.png",type=c("cairo"), height = 10, width = 11, units = 'in', res = 400)
Heavy2 %>%
  ggplot(aes(UMAP_1, UMAP_2))+
  geom_point(shape=16, color="darkgray", size=0.5, alpha=0.5)+
  geom_point(data=Heavy2 %>% filter(clonal_abundance>2), shape=19, alpha=0.5, color="black", size=1)+
  geom_point(data=Heavy2 %>% filter(clonal_abundance>2), shape=16, size=1, aes(color=total_mut))+
  geom_line(data=Heavy2 %>% filter(clone_per_replicate=="108-N00044-Ja001"), size=0.2)+
  scale_color_distiller(palette = 'RdYlBu')+
  #scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 11, name = "Spectral")))(50))+
  facet_wrap(~Source)
dev.off()

Heavy2 %>%select(clone_per_replicate, clonal_abundance, lv_dist) %>% arrange(desc(clonal_abundance)) %>% filter(clonal_abundance <700)
#Plot of proposed activated naive B cells defined by mutations >2 and naive annotations on GEX.
png("/well/combat/users/vkh192/repertoire/out/umap_act_naive.png",type=c("cairo"), height = 6, width = 8, units = 'in', res = 400)
Heavy %>%
  ggplot(aes(UMAP_1, UMAP_2, color=(total_mut>2 & grepl("naive", name))))+
  geom_point(alpha=0.7, size=0.7)+
  scale_color_manual(values=c("gray", "red"))
dev.off()

#Facet activated naive plots by c_gene
png("/well/combat/users/vkh192/repertoire/out/umap_act_naive_facet.png",type=c("cairo"), height = 6, width = 10, units = 'in', res = 400)
Heavy %>%
  mutate(c_gene=gsub("\\*.*", "", c_gene)) %>%
  ggplot(aes(UMAP_1, UMAP_2, color=(total_mut>2 & grepl("naive", name))))+
  geom_point(alpha=0.7, size=0.4)+
  scale_color_manual(values=c("gray", "red"))+
  facet_grid(~c_gene)
dev.off()

