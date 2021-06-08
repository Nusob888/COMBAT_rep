## This code is adapted from BCellR (unpublished R package for single cell analysis of B cell VDJ immune profiling data). 
## This code will generate clonotypes according to the following procedure: 
## CD-HIT clustering of VJ gene-pair sequences followed by within-cluster local minimum thresholding of pairwise CDR3 normalised levenstein distances.
## Where sample depth is low, automatic thresholding will not be robust, therefore a heuristic threshold should be set by manual inspection of sequences where possible. 

## Author: Bo Sun
## Group: Bashford-Rogers
## Usage: If code is to be reproduced/adapted for other works, please aknowledge code authorship, or cite the BCellR package (if package has been formally released at the time of reproduction) 

#'@importFrom Biostrings DNAStringSet
#'@importFrom Biostrings writeXStringSet
#'@importFrom Biostrings AAStringSet
#'@importFrom igraph graph.data.frame
#' @importFrom magrittr %>%
#' @importFrom  tidyr separate
#' @importFrom  tidyr gather
#' @import dplyr
#' @import  parallel 
#' @import data.table
#' @import stringdist

#cdhit function, calling on cdhit functions in src
cdhit <- function (seqs, options, name = "CD-Hit", showProgress = interactive()) 
{
  seqs <- Biostrings::DNAStringSet(seqs)
  options$i <- tempfile()
  writeXStringSet(seqs, options$i)
  on.exit(unlink(options$i))
  switch(class(seqs), AAStringSet = cdhitC(options, name, showProgress) + 
           1, DNAStringSet = BcellR:::cdhitestC(options, name, showProgress) + 
           1, stop("seqs must be either AAStringSet or DNAStringSet"))
}

#function to calculate levenshtein distances of cdr3s after VJ clustering
stringdistclusters <- function(clusterlist){
  numCores <- detectCores()
  clusterlist <- lapply(clusterlist, function(x){
    x$length <- nchar(x$junction_nt)
    dist <- stringdist::stringdistmatrix(x$junction_nt, x$junction_nt, method="lv")
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
    x <- dist %>% select(-s1_length, -s2_length)
  })
  edgelist <- data.table::rbindlist(clusterlist) 
  return(edgelist)
}

#function to calculate levenshtein distances of cdr3s after VJ clustering
igraphclusters <- function(cdr3dist){
  numCores <- detectCores()
  #create list per replicate
  cdr3list <- split(cdr3dist, cdr3dist$COMBAT_ID_Time)

  cdr3list <- mclapply(cdr3list, function(x){
    distgraph <- igraph::graph.data.frame(x[,1:2])
    links <- data.table::data.table(sequence_number=unique(unlist(x[,1:2])), clone_per_replicate= paste(igraph::clusters(distgraph)$membership, unique(unlist(x[,c("COMBAT_ID_Time")])), sep = "-"))
    x<- x%>% left_join(links, by= c("sequence_number"))
  }, mc.cores = numCores) 

  clones <- data.table::rbindlist(cdr3list) 
  
  #create list per baseID
  clonelist <- split(clones, clones$COMBAT_ID)

  clonelist <- mclapply(clonelist, function(x){
    distgraph <- igraph::graph.data.frame(x[,1:2])
    links <- data.table::data.table(sequence_number=unique(unlist(x[,1:2])), clone_per_baseID= paste(igraph::clusters(distgraph)$membership, unique(unlist(x[,c("COMBAT_ID")])), sep = "-"))
    x<- x%>% left_join(links, by= c("sequence_number"))
  }, mc.cores = numCores) 
  
  clones <- data.table::rbindlist(clonelist) 
  
  #create public clone list
  
    distgraph <- igraph::graph.data.frame(clones[,1:2])
    links <- data.table::data.table(sequence_number=unique(unlist(clones[,1:2])), clone_global= paste(igraph::clusters(distgraph)$membership, "ALLCOMBAT", sep="-"))
    clones <- clones %>% left_join(links, by= c("sequence_number"))
  
  return(clones)
}

#' @export
#Function to remove last nucleotide of full length sequences not divisible by three
correct10Xfullseq <- function(x){
  x<- x %>% 
    dplyr::mutate(sequence = ifelse((!as.numeric(nchar(sequence)) %% 3 == 0), str_sub(sequence, end=-2), sequence))
  return(x)
}

#' @export
#function to get cdr3 lv distances after 80% cdhit clustering of VJ
Get10XclonesCOMBAT <- function(x) {

  x <- Heavy %>% 
    group_by(sample, sequence) %>% 
    add_tally(name="clonesize") %>% 
    ungroup()
  
  #remove R02 values from longitudinal samples
  #concatenate V and J segments into a new vector
  vj <- x %>%  
    tidyr::unite(v_j_pair, c("v_region", "j_region"), sep="", remove=FALSE)
  vjseq <- Biostrings::DNAStringSet(vj$v_j_pair)
  names(vjseq) <- as.numeric(vj$sequence_number)
  
  #Set CDhit opts, for now not changeable to user directly. Can incorporate if needed
  cdhitOpts <- setNames(as.list(c(0.60, 0, 40)), c("c", "M", "AL"))
  cdhitOpts <- lapply(cdhitOpts, as.character)
  
  #run cdhit to cluster vj segments of 80% identity
  cdhitclusters<- cdhit(vjseq, cdhitOpts, 'Preclustering')
  #extract corresponding whole vdj
  vdj <- data.table::as.data.table(cbind(names(vjseq), unlist(cdhitclusters),x$sequence, x$junction_nt))
  colnames(vdj) <- c("sequence_number", "cdhit_cluster", "sequence", "junction_nt")
  vdj$sequence_number <- as.numeric(vdj$sequence_number)
  
  #create list of clusters
  clusterlist <- split(vdj, vdj$cdhit_cluster)
  
  cdr3dist<- stringdistclusters(clusterlist)
  head(cdr3dist)  


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
  cdr3dist$sequence_number <- as.character(cdr3dist$sequence_number)
  cdr3dist$seq_to <- as.character(cdr3dist$seq_to)
  x$sequence_number <- as.character(x$sequence_number)
  
  cdr3dist <- cdr3dist %>% left_join(x, by = "sequence_number") 
  
  #add clones
  cdr3dist<- igraphclusters(cdr3dist)
head(cdr3dist)  

cdr3dist %>% select(sample)
  edgelist <- cdr3dist[,c("sequence_number", "seq_to", "sequence")]
  nam <- cdr3dist %>% ungroup() %>% group_by(sequence_number) %>% dplyr::slice(1) %>% ungroup()
  
  #Get clone abunances
  nam2 <- nam %>% 
    dplyr::group_by(barcode_seq,COMBAT_ID_Time) %>%
    slice(1) %>% 
    dplyr::ungroup() %>%
    dplyr::group_by(sample) %>% 
    dplyr::add_tally(name= "tot_samp") %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(clone_per_replicate, COMBAT_ID_Time) %>% 
    dplyr::add_tally(name = "clonal_abundance") %>% 
    dplyr::ungroup() %>%
    dplyr::select(barcode_seq, clonal_abundance)
  
  nam <- nam %>% left_join(nam2, by = "barcode_seq")
  
  out<- list(nam, edgelist)
  names(out) <- c("Heavy", "edgelist")

  return(out)
}
















