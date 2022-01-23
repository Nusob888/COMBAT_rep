##create combat meta
#scp -r /Users/bosun/Documents/GitHub/COMBAT_rep/08_CZI_table_deposit/barcodes_in_analysis.txt  vkh192@rescomp.well.ox.ac.uk:/well/combat/users/vkh192/repertoire/
#scp -r /Users/bosun/Documents/GitHub/COMBAT_rep/08_CZI_table_deposit/BCR_barcodes.txt vkh192@rescomp.well.ox.ac.uk:/well/combat/users/vkh192/repertoire/

library(tidyverse)

BCR <- data.table::fread("/well/combat/datasets/CBD-10X-00009/bcr_release_v001_chainmerge_table.csv.gz")
BCR_barcodes <- data.table::fread("/well/combat/users/vkh192/repertoire/BCR_barcodes.txt")
  
TCR <- data.table::fread("/well/combat/datasets/CBD-10X-00008/tcr_chain_information.tsv.gz")
TCR_barcodes <- data.table::fread("/well/combat/users/vkh192/repertoire/barcodes_in_analysis.txt")
head(BCR)
head(TCR)
nrow(TCR)

BCRtrim <- BCR %>%
  select(barcode_seq, contig_qc_HC, umis_HC, 
         functionality_HC, v_call_HC,v_score_HC, 
         j_call_HC, j_score_HC, junction_aa_HC, total_mut_HC, s_mut_HC, r_mut_HC,
         c_gene_HC, clone_per_replicate_HC, clone_global_HC, clonal_abundance_HC, 
         locus_LC, contig_qc_LC, umis_LC, functionality_LC, v_call_LC, v_score_LC, 
         j_call_LC, j_score_LC, junction_aa_LC, total_mut_LC, s_mut_LC, r_mut_LC, 
         c_gene_LC)
BCRtrim <- BCRtrim %>%
  rename(barcode_id=barcode_seq)
colnames(TCR)

TCRtrim <- TCR %>%
  select(barcode_id, chain_composition, clone_ID, clone_count, clone_proportion, 
         contains_unproductive, doublet, chain_TRA, v_gene_TRA, d_gene_TRA, j_gene_TRA, 
         c_gene_TRA, productive_TRA, cdr3_TRA, umis_TRA, chain_TRA2, chain_TRA2, v_gene_TRA2,
         d_gene_TRA2, j_gene_TRA2, c_gene_TRA2, productive_TRA2, cdr3_TRA2, umis_TRA2,
         chain_TRB, v_gene_TRB, d_gene_TRB, j_gene_TRB, c_gene_TRB, productive_TRB, 
         cdr3_TRB2, umis_TRB2, chain_TRB2, v_gene_TRB2, d_gene_TRB2, j_gene_TRB2, c_gene_TRB2, productive_TRB2, 
         cdr3_TRB2, umis_TRB2)

BCRexample <- BCRtrim[1:20, ]
TCRexample <- TCRtrim[1:20, ]

write.csv(BCRexample, "/well/combat/users/vkh192/repertoire/BCR_example.csv")
write.csv(TCRexample, "/well/combat/users/vkh192/repertoire/TCR_example.csv")

write.csv(BCRtrim, "/well/combat/users/vkh192/repertoire/BCR_VDJ_CZI.csv")

write.csv(TCRtrim, "/well/combat/users/vkh192/repertoire/TCR_VDJ_CZI.csv")


BCRtrim <- data.table::fread("/well/combat/users/vkh192/repertoire/BCR_VDJ_CZI.csv") %>% select(-V1)

TCRtrim <- data.table::fread("/well/combat/users/vkh192/repertoire/TCR_VDJ_CZI.csv")%>% select(-V1)

BCRtrim$used_in_analysis <- BCRtrim$barcode_id %in% BCR_barcodes$x
TCRtrim$used_in_analysis <- TCRtrim$barcode_id %in% TCR_barcodes$x

summary(TCRtrim$used_in_analysis)
length(TCR_barcodes$x)

write.csv(BCRtrim, "/well/combat/projects/repertoire/CZI_deposit/BCR_VDJ_CZI.csv")

write.csv(TCRtrim, "/well/combat/projects/repertoire/CZI_deposit/TCR_VDJ_CZI.csv")
