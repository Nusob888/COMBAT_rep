## This code is adapted from BCellR (unpublished R package for single cell analysis of B cell VDJ immune profiling data). 

## This code will assign contig QC metrics for downstream parsing based on the following procedure:
## All contigs are assigned a rank order per cell barcode. Where there are more than 1 contigs per Heavy or Light chain locus, the top two ranked contigs are selected and a rank UMI ratio is calculated (rank2/rank1). 
## The knee point of the log ranked UMI ratios is then taken as the threshold for assigning contig confidence as a VDJ singlet vs doublet. 
## Where the clonal identity (identical VJpairing and CDR3 nucleotide sequence) is matched for rank1 and rank2 contigs that are IGHM/IGHD, both contigs are retained. 
## Whilst double light chains are potentially physiologically relevant in 1% of B cells, we have opted to keep only the chains that pass QC. This is due to the inability to reliably discriminate between ambient RNA light chain contamination from true light chain expression.

## Author: Bo Sun
## Date: 2020-08-11
## Group: Bashford-Rogers
## Usage: Please do not reproduce this code outside of the COMBAT environment until official publication. 
## If code is to be reproduced/adapted for other works, please aknowledge code authorship, or cite the BCellR package (if package has been formally released at the time of reproduction) 

#' @import stringdist
#' @import data.table
#' @import parallel


#' @export
Parsechains <- function(data, filter_productive=TRUE, format="IMGT", sampleid = "gplex"){
  numCores <- detectCores()-1
  #Create locus column for both IMGT and 10X formats
  #This is required to bypass "Multi chain" assignment by 10X which doesn't allow group splitting
  if (missing(sampleid)){
    stop("please specify the sample pool annotation column")
  } else if (missing(filter_productive)){
    stop("please specify if you would like to filter for productive sequences")
  } else if (filter_productive==TRUE){
    if(missing(format)){
      stop("oops, you didn't specify the input format: IMGT? or 10X?")
    }else if (format=="IMGT"){
      data <- data %>% 
        mutate(barcode_seq = paste0(barcode, "-",get(sampleid), sep="")) %>%
        filter(grepl("^productive", `functionality`)) %>%
        dplyr::select(barcode_seq, v_call, j_call, junction_aa, locus, umis, contig_id, c_gene) %>% #new
        group_by(barcode_seq, locus) %>%
        mutate(umi_ranks = order(order(umis,contig_id, decreasing=TRUE))) %>% #add umi ranks per contig
        ungroup() %>% 
        filter(umi_ranks <=2) #Here we take the top 2 ranked contigs. I need to extend the function to handle more than 2 contigs in case of clonally related contigs
    }else if (format=="10X"){
      df <- data %>% 
        mutate(barcode_seq = paste0(barcode, "-",get(sampleid), sep="")) %>%
        filter(grepl("^True|^TRUE", productive)) %>%
        mutate(locus = ifelse(grepl("IGH", v_gene), "Heavy", ifelse(grepl("IG[KL]", v_gene), "Light", "NA"))) %>%        
        dplyr::select(barcode_seq, v_gene, j_gene, cdr3, locus, umis, contig_id, c_gene) %>% #new
        group_by(barcode_seq, locus) %>%
        mutate(umi_ranks = order(order(umis,contig_id, decreasing=TRUE))) %>% #add umi ranks per contig
        ungroup() %>% 
        filter(umi_ranks <=2)
    }
  }else if (filter_productive == FALSE){
    data
  }
  
  #Add VDJ_identities
  if (format == "IMGT"){
    data$VDJ_group <- data %>% 
      group_indices(barcode_seq, v_call, j_call, junction_aa)
  } else if (format == "10X"){
    data$VDJ_group <- data %>% 
      group_indices(barcode_seq, v_gene, j_gene, cdr3) 
  }
  
  #Get ratio cutoff using mixture model fit
  ranks <- data %>% 
    select(barcode_seq,contig_id, umis , locus, umi_ranks) %>% 
    group_by(barcode_seq, locus) %>%
    add_tally(name="n_contigs") %>% 
    filter(n_contigs>=2) %>% 
    summarise(rank_ratio = umis[umi_ranks == 2] / umis[umi_ranks == 1]) %>%
    ungroup() %>%
    mutate(contig_ranks = row_number(desc(rank_ratio))) %>% 
    arrange(desc(contig_ranks))
  
  # find inflection point - code taken from Kallisto package 
  get_inflection <- function(df, lower = 0.1) {
    log_total <- log_rank <- rank_ratio <-  NULL
    df_fit <- df %>% 
      dplyr::filter(rank_ratio > lower) %>% 
      transmute(log_total = log10(rank_ratio),
                log_rank = log10(contig_ranks))
    d1n <- diff(df_fit$log_total)/diff(df_fit$log_rank)
    right.edge <- which.min(d1n)
    10^(df_fit$log_total[right.edge])
  }
  cutoff <<- get_inflection(ranks)
  
  #split by barcodes suffixed by sample. This ensures no clashing barcodes
  split <-  data %>% 
    select(barcode_seq, locus, umis, umi_ranks, contig_id, VDJ_group, c_gene) %>% 
    group_by(barcode_seq) %>% 
    group_split()
  
  #Fun filtercontigsscript
  split <- mclapply(split, function(x){filtercontigs(x, cutoff)},mc.cores=numCores)
  split <- data.table::rbindlist(split, use.names = TRUE)
  split <- split %>% filter(!is.na(barcode_seq)) %>% select(-c_gene, VDJ_group, umi_ranks)
  return(split)
}


#' @export
filtercontigs<- function(test, cutoff){
  
  #Parse chains
  if(any(grepl("IGH|Heavy", test$locus))){
    Heavy <- test %>% filter(grepl("IGH|Heavy", locus))
  }
  
  if(any(grepl("IGK|IGL|Light", test$locus))){
    Light <- test %>% filter(grepl("IGK|IGL|Light", locus))
  }
  
  #Process Heavy chain
  if(exists("Heavy")){
    if(!any(Heavy$umi_ranks== 2)){
      heavy <- Heavy[c(Heavy$umi_ranks == 1),]#select rank1 if rank2 ratio less than < 0.5 (rank1 umis is double of rank2)
      heavy$contig_qc <- "singleton"
    }else{ 
      for(i in 1:nrow(Heavy)){
        if (Heavy$umi_ranks[i] == 1){#----
          if(Heavy$c_gene[i] == "IGHM|IGHD"){#conditional selection of rank1 IGHMs
            if(Heavy[c(Heavy$umi_ranks== 1), c("VDJ_group")] == Heavy[c(Heavy$umi_ranks== 2), c("VDJ_group")]){
              heavy <- Heavy %>% 
                filter(grepl(paste(Heavy$VDJ_group[i]), VDJ_group)) %>% #filter on VJ and cdr3 matches
                filter(grepl("IGHM|IGHD", c_gene)) #filter VJ/cdr3 matches that are IGHM or IGHD, thereby for top IGHMs, keeping paired IGHD expression
              heavy$contig_qc <- "passed_qc"
            }else if((Heavy[c(Heavy$umi_ranks == 2),c("umis")]/Heavy[c(Heavy$umi_ranks == 1), c("umis")]) <cutoff){
              heavy <- Heavy %>% 
                filter(grepl(paste(Heavy$VDJ_group[i]), VDJ_group)) %>% #filter on VJ and cdr3 matches
                filter(grepl("IGHM|IGHD", c_gene)) #filter VJ/cdr3 matches that are IGHM or IGHD, thereby for top IGHMs, keeping paired IGHD expression
              heavy$contig_qc <- "passed_qc"
            }else{
              heavy <- Heavy
              heavy$contig_qc <- "potential_doublet"
            }
          }else if(test$c_gene[i] != "IGHM"){# conditional selection of non_IGHM rank1s
            if(Heavy[c(Heavy$umi_ranks== 1), c("VDJ_group")] == Heavy[c(Heavy$umi_ranks== 2), c("VDJ_group")]){
              heavy <- Heavy[c(Heavy$umi_ranks == 1),]
              heavy$contig_qc <- "passed_qc"
            }else if((Heavy[c(Heavy$umi_ranks == 2),c("umis")]/Heavy[c(Heavy$umi_ranks == 1), c("umis")]) <cutoff){
              heavy <- Heavy[c(Heavy$umi_ranks == 1),]
              heavy$contig_qc <- "passed_qc"
            }else{
              heavy <- Heavy
              heavy$contig_qc <- "potential_doublet"
            }
          }
        }
      }
    }
  }
  
  #Process light chain
  if(exists("Light")){
    if(!any(Light$umi_ranks== 2)){
      light <- Light[c(Light$umi_ranks == 1),] #select rank1 if rank2 ratio less than < 0.5 (rank1 umis is double of rank2)
      light$contig_qc <- "singleton"
    }else{
      for(i in 1:nrow(Light)){
        if (Light$umi_ranks[i] == 1){
          if((Light[c(Light$umi_ranks == 2), c("umis")] / Light[c(Light$umi_ranks == 1), c("umis")]) < cutoff){
            light <- Light[c(Light$umi_ranks == "1"),]
            light$contig_qc <- "passed_qc"
          }else{
            light <- Light #need to add check for alternative lambda or kappa
            light$contig_qc <- "potential_doublet"
            
          }
        } 
      }
    }
  }
  
  #Sort objects
  if(exists("heavy") & exists("light")){
    merge <- mapply(c, heavy, light, SIMPLIFY=FALSE)
    return(merge)
  }else if(exists("heavy") & isFALSE(exists("light"))){
    return(heavy)
  }else if(exists("light") & isFALSE(exists("heavy"))){
    return(light)
  }
}
