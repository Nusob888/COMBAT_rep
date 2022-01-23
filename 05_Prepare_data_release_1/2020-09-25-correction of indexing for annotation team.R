#Script to check annotation for annotation team
##Author: Bo Sun
##Lab: Bashford Rogers
##Institute: University of Oxford

#Load annotation team data
anno <- readr::read_csv("/well/combat/users/vkh192/repertoire/data/b_cell_reannotation_all_barcodes.csv")
banno <- anno %>% filter(subset=="D")
##Load BCR rds
clones <- readRDS("/well/combat/users/vkh192/repertoire/data/new_clones.rds")
Heavy <- clones$Heavy

#load new annotations previously made by me
newanno <- readr::read_csv("/well/combat/users/vkh192/repertoire/data/newanno.csv") %>% select(-name)
contig_qc <- readr::read_csv("/well/combat/projects/repertoire/scbcr_table-001/contig_qc-001.csv")

#join datasets
Heavy <- Heavy %>% left_join(newanno, by=c("barcode_seq"="X1"))
Heavy <- Heavy %>% select(barcode_seq, COMBAT_ID_Time, new_clusters)
Heavy <- Heavy %>% left_join(scbcr %>% select(-X1) %>% group_by(barcode_seq) %>% slice(1), by=c("barcode_seq", "COMBAT_ID_Time"))
nrow(Heavy)
Heavy <- merge(anno, Heavy, by.x="barcode", by.y="barcode_seq", all.x=TRUE)
head(Heavy)
nrow(Heavy)
head(as.data.frame(Heavy)) #note here that the parse_chains.y from the annotation team merge seems to have indexed incorrectly

Heavy[!(anno$barcode %in% Heavy$barcode_seq),] %>% filter(subset=="D") %>% nrow

#Create correctly indexed file for annotation team

#Select names as per annotation team file with addition of my putative annotations
Heavy <- Heavy %>% select(barcode, cluster_id, cellID, subset, name, new_clusters, pseudobulk, major_cell_type, ADT_cluster_no, ADT_cluster_name, parse_chains.y, repertoire_naive_activated)
Heavy<- Heavy %>% rename(parse_chains=parse_chains.y)

#Sanity check indexing
head(as.data.frame(Heavy)) #looks good

#Select on B cells (B) and sanity check for missing annotations
nrow(banno) #57438 B cells in subset D
nrow(clones$Heavy) #49516 B cells present in VDJ data
57438-49516 #7922 missing from VDJ

Heavy %>% filter(subset == "D") %>% filter(is.na(new_clusters)) %>% nrow #7922 missing annotations in my new assignments 
Heavy %>% filter(subset == "D") %>% filter(is.na(name)) %>% nrow #No missing annotations in the old annotations
Heavy %>% filter(subset == "D") %>% filter(is.na(parse_chains)) %>% nrow #7922 missing values for chain QC check

#Save correctly indexed file
write.csv(Heavy, "/well/combat/users/vkh192/repertoire/data/b_cell_reannotation_all_barcodes_boedit.csv")
nrow(anno)
nrow(Heavy)






##Sofija data
caspr2db <- data.table::fread("data/2020-07-26-Uptodate_CASPR2db_noseq-v2.txt", na.strings = c("NA", ""), stringsAsFactors = FALSE)


caspr2db %>%
  mutate(mAb_ID = Sequence.ID) %>%
  mutate(mAb_ID= factor(mAb_ID, levels=c("E04","E04R", "E06", "E07", "E07R", "H01", "H01R", "H02" ,"H02R",
                                         "E08", "E08R", "E09", "E09R"))) %>%
  filter(!grepl("R", mAb_ID))%>%
  select(mAb_ID,Subtype, Sofija_human_end_titer_1, Sofija_mouse_end_titer_1,`Specificity_Bo_C2/C4_swap`) %>%
  reshape2::melt(id.vars=c("mAb_ID", "Subtype", "Specificity_Bo_C2/C4_swap")) %>%
  unique %>%
  filter(!is.na(value)) %>%
  mutate(variable = ifelse(grepl("human", variable), "Human", "Mouse")) %>%
  ggplot(aes(variable, value, group=mAb_ID))+
  geom_line(aes(group=mAb_ID), alpha=0.5, linetype="dashed")+
  geom_point(shape=23,size=2, color="black",aes(fill=mAb_ID), fill="black")+
  scale_y_log10()+
  facet_grid( vars(`Specificity_Bo_C2/C4_swap`),vars(mAb_ID))+
  ylab("log10 end-point titer (1:X)")+
  annotation_logticks(sides = "l", size=0.1, outside = FALSE) +
  xlab("")
colnames(caspr2db)

Disc <- ggmsa("data/_out.210806011828332eNrjl3GnEL9g1owqYhTcClsfnormal.fasta", 35, 185, color = "Clustal", font = NULL, seq_name = T ) + geom_msaBar()
Disc
Lam3 <- ggmsa("data/_out.210806011828332eNrjl3GnEL9g1owqYhTcClsfnormal.fasta", 799, 963, color = "Clustal", font = NULL, char_width = 0.5, seq_name = T )+ geom_msaBar()
Lam3

