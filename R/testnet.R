#sanity checking edgelist
test <- clones
test$edgelist %>% left_join(test$Heavy, by = c("sequence_number")) %>% arrange(sequence_number)

length(c(test$edgelist$sequence_number, test$edgelist$seq_to) %>% unique) 

summary((c(test$edgelist$sequence_number, test$edgelist$seq_to) %>% unique) %in% ((test$Heavy$sequence_number) %>% unique))

length((test$Heavy$sequence_number %>% unique))

pdf("/well/combat/users/vkh192/repertoire/out/cdr3dist.pdf", height = 20, width = 10)
ggplot(cdr3dist, aes(x=norm_lvdist))+
  geom_histogram(binwidth=0.01) +
  geom_vline(xintercept = cutoff)
dev.off()

library(GGally)
library(ggnetwork)
library(network)
library(sna)
test <- clones

test$edgelist <- merge(test$edgelist, test$Heavy, c("sequence_number", "seq_to"), all.x=TRUE)
test$distinct <- test$edgelist %>% group_by(sequence_number) %>% dplyr::slice(1) %>% ungroup %>% arrange(sequence_number)

test$edgelist %>% ungroup() %>% select(clone, junction_nt, v_call, j_call, total_mut) %>% unique

#sanity check
length(test$distinct$sequence_number)
test$edgelist %>% ungroup() %>% group_by(sequence_number, seq_to) %>% unique() %>% nrow

set.seed(35)
n = network::network(test$edgelist[, 1:2 ], directed = TRUE) # directed graph




n %v% "degree" = sna::degree(n)
n %v% "Patient" = test$distinct$sample
n %v% "vfam" = test$distinct$vfam
n %v% "clonesize" = test$distinct$clonesize
n %v% "mutations" = test$distinct$total_mut
n %v% "clone" = test$distinct$clone
n %v% "Source" = test$distinct$source
n %v% "c_gene" = gsub("\\*.*","", test$distinct$c_gene)
n %v% "barcode" = test$distinct$barcode
n %v% "umis" = test$distinct$umis


set.seed(1)
net <- ggplot(ggnetwork::fortify(n, arrow.gap = 0), aes(x, y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(fill = Patient, size=mutations), shape = 21) +
  theme_blank() 
pdf("/well/combat/users/vkh192/repertoire/out/prelimbnet_release1.pdf", height = 20, width = 20)
net
dev.off()



pdf("/well/combat/users/vkh192/repertoire/out/net_facet_source.pdf", height = 20, width = 20)
ggplot(ggnetwork::fortify(n, arrow.gap = 0), aes(x, y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(fill = as.numeric(mutations), size=clonesize), shape = 21) +
  theme_blank() +
  facet_wrap(~Source)
dev.off()

pdf("/well/combat/users/vkh192/repertoire/out/net_facet_source_iso.pdf", height = 20, width = 20)
ggplot(ggnetwork::fortify(n, arrow.gap = 0), aes(x, y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(fill = c_gene, size=degree), shape = 21) +
  theme_blank() +
  facet_wrap(~Source)
dev.off()

pdf("/well/combat/users/vkh192/repertoire/out/net_facet_source_iso_PC.pdf", height = 20, width = 20)
ggplot(ggnetwork::fortify(n, arrow.gap = 0), aes(x, y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(fill = c_gene, size=umis), shape = 21) +
  theme_blank() +
  facet_wrap(~Source)
dev.off()






scp -r vkh192@rescomp.well.ox.ac.uk:/well/combat/users/vkh192/repertoire/out/prelimbnet_release1.pdf /Users/bosun/Documents/DPhil/Projects/COMBAT/Repertoire/output

network::list.vertex.attributes(n)
vertexnames<- network::network.vertex.names(n)
vertexpatient<- network::get.vertex.attribute(n, "Patient", unlist = TRUE)
degree<- network::get.vertex.attribute(n, "degree", unlist = TRUE)
clone<- network::get.vertex.attribute(n, "clone", unlist = TRUE)
barcode<- network::get.vertex.attribute(n, "barcode", unlist = TRUE)
clonesize<- network::get.vertex.attribute(n, "clonesize", unlist = TRUE)

sanity<- data.table::data.table(cbind(vertexnames, degree, vertexpatient, clone, barcode, clonesize))
sanity$degree <- as.numeric(sanity$degree)
head(sanity %>% dplyr::arrange(desc(clonesize)))

head(clones$Heavy %>% ungroup %>% filter(clone == "27-COMBAT0070R01") %>%select(sequence_number, sample, clone), 30)
head(sanity %>% ungroup %>% filter(clone == "27-COMBAT0070R01") %>%select(vertexnames, vertexpatient, clone), 30)

head(test$Heavy %>% filter(sequence_number == "111372") %>% select(clone, COMBATID))
