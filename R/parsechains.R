#wrapper for Filter_VDJ_contig_HL.R

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
        mutate(barcode_by_sample = paste0(barcode, get(sampleid), sep="")) %>%
        filter(grepl("^productive", `V-DOMAIN Functionality`)) %>%
        mutate(locus = gsub("Homsap ", "", `V-GENE and allele`)) %>%
        mutate(locus = gsub("V.*", "", locus)) %>%
        group_by(barcode_by_sample, locus) %>%
        mutate(umi_ranks = order(order(umis,contig_id, decreasing=TRUE))) %>% #add umi ranks per contig
        ungroup() %>% 
        filter(umi_ranks <=2) #Here we take the top 2 ranked contigs. I need to extend the function to handle more than 2 contigs in case of clonally related contigs
    }else if (format=="10X"){
      data <- data %>% 
        mutate(barcode_by_sample = paste0(barcode, get(sampleid), sep="")) %>%
        filter(grepl("^productive", v_gene)) %>%
        mutate(locus = gsub("V.*", "", v_gene)) %>%
        group_by(barcode_by_sample, IMGTchain) %>%
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
      group_indices(barcode_by_sample,`V-GENE and allele`, `J-GENE and allele`, `AA JUNCTION`)
  } else if (format == "10X"){
    data$VDJ_group <- data %>% 
      group_indices(barcode_by_sample, v_gene, j_gene, cdr3) 
  }
  
  #split by barcodes suffixed by sample. This ensures no clashing barcodes
  split <-  data %>% 
    dplyr::select(barcode_by_sample, locus, umis, umi_ranks, contig_id, VDJ_group, c_gene) %>% 
    group_by(barcode_by_sample, locus) %>% 
    group_split()
  
  #Fun filtercontigsscript
  split <- mclapply(split, function(x){filtercontigs(x)},mc.cores=numCores)
  split <- data.table::rbindlist(split, use.names = TRUE)
  
  return(split)
}

