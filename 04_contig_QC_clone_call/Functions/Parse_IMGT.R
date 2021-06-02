#Code to parse IMGT outputs for BCR

#' @importFrom magrittr %>%
#' @importFrom  tidyr separate
#' @importFrom  tidyr gather
#' @import dplyr
#' @import  parallel 
#' @import data.table
#' @import stringr
#' @import Peptides


#Function taken from bcRep credit: Julia.Bischof@uksh.de
readIMGT<-function(data,filterNoResults=TRUE){
  if(substr(strsplit(data,split="/")[[1]][length(strsplit(data,split="/")[[1]])],1,2) %in% c("3_","6_")){
    newcolnames<-gsub(" |-","_",strsplit(readLines(data)[1],split="\t")[[1]])
    newcolnames<-gsub("'","",newcolnames)
    
    if(length(strsplit(readLines(data)[2],split="\t")[[1]])==3){
      temp<-lapply(readLines(data),function(x){strsplit(x,split="\t")[[1]]})
      for(i in 1:length(temp)){
        if(length(temp[[i]])==3){
          temp[[i]]<-c(temp[[i]],rep(" ",(length(newcolnames)-3)))
        }
      }
      tab<-t(data.frame(temp[2:length(temp)], stringsAsFactors = F))
      rownames(tab)<-as.numeric(seq(1, (length(temp)-1), 1))
    }else{
      tab<-data.frame(matrix(NA, nrow=0,ncol=length(newcolnames)))
      tab<-rbind(tab,read.table(data,fill=TRUE,header=F,skip=1,sep="\t",as.is=TRUE,na.strings=c("","NA")))
    }
    colnames(tab)<-newcolnames
  }else{
    tab<-read.table(data,fill=TRUE,header=TRUE,sep="\t",as.is=TRUE,check.names=F, na.strings=c("","NA"))
    colnames(tab)<-gsub(" |-","_",gsub("'","",colnames(tab)))
  }
  tab<-tab[,-ncol(tab)]
  if(filterNoResults==T){
    out<-grep("No results|No rearrangement found|unproductive",tab[,3])
    if(length(out)>0){
      tab<-tab[-out,]
      cat("...", length(out), "sequences were excluded\n")
    }else{
      cat("... no sequences were excluded\n")
    }      
  }
  return(tab)
}

#Function to parse hotspot txt
Get_hotspot_df <- function(hotspot_tab) {
  hotspots <- hotspot_tab %>%
    select(-c(1:5, 11)) %>%
    group_by(V_GENE_and_allele) %>%
    unique()
  colnames(hotspots) <- c("V_GENE_and_allele", "WA", "TW", "RGYW", "WRCY")
  hotspots <- hotspots %>% tidyr::gather(Hotspot, Positions, "WA", "TW", "RGYW", "WRCY")
  hotspots$Positions <- strsplit(hotspots$Positions, '\\|')
  for (i in 1:length(hotspots$Positions)) {
    hotspots$Positions[[i]]<- reshape2::colsplit(hotspots$Positions[[i]], "\\,", names=c("Motif", "Position")) %>%
      tidyr::separate(Position, c("Position", "Region"), sep = "\\(") %>%
      dplyr::mutate(Region = gsub(")", "", Region))
  }
  return(hotspots)
}


#Function to parse a list of IMGT mutation objects
parsemut<- function(x){
  blank <- setNames(data.table::data.table(matrix(ncol = 7, nrow = 1)), c("nt_change", "nt_position", "aa_change","aa_similarity", "aa_position", "codon_from", "codon_to"))
  if (anyNA(x)){
    x <- blank
  }else{
    x <- data.table::setDT(data.table::tstrsplit(x, ",", names=c("nt_change", "Position"))) %>%
      .[, c("nt_position") := stringr::str_extract(nt_change, "[[:digit:]]+")] %>%
      .[, nt_change := stringr::str_remove(nt_change, "[[:digit:]]+")] %>%
      .[, c("aa_change", "Position") := data.table::tstrsplit(Position, ";")] %>%
      .[, c("aa_change", "aa_similarity") := data.table::tstrsplit(aa_change, "\\(|\\)")] %>%
      .[, c("aa_position") := stringr::str_extract(aa_change, "[[:digit:]]+")] %>%
      .[, aa_change := stringr::str_remove(aa_change, "[[:digit:]]+")] %>%
      .[, Position := stringr::str_remove(Position, "\\s*\\[[^\\)]+\\]")] %>%
      .[, c("blank", "aa_from", "codon_from", "position", "codon_to") := data.table::tstrsplit(Position, "[[:space:]]")] %>%
      .[, c("nt_change", "nt_position", "aa_change","aa_similarity", "aa_position", "codon_from", "codon_to")] %>%
      .[, aa_change := ifelse(!grepl(">|NA", aa_change), "Silent", aa_change)]
  }
  return(x)
}


ParseIMGTmut <- function(test){
  test <- data.table::as.data.table(test)
  numCores <- parallel::detectCores()
  cols <- c("Sequence_number", "FR1_IMGT", "CDR1_IMGT", "FR2_IMGT", "CDR2_IMGT", "FR3_IMGT") #we do not analyse CDR3 region here due to unreliable alignment
  subset <- test[, ..cols]
  for (i in 2:length(subset)){
    subset[[i]] <- strsplit(subset[[i]], '\\|') #Creates column of lists of observed mutations
  }
  for (i in 2:length(subset)){
    subset[[i]] <- parallel::mclapply(subset[[i]], parsemut, mc.cores=numCores) #Parses each list item with parsemut2
  }
  return(subset)  
}

compileIMGTdfs <- function(imgtdfs){
  cat("Parsing IMGT output")
  IMGTsum<- imgtdfs[["1_Summary.txt"]] %>%
    dplyr::select("Sequence_number","barcode", "contig", "sample", "V_DOMAIN_Functionality", "V_GENE_and_allele", "V_REGION_score", "V_REGION_identity_%", 
                  "V_REGION_potential_ins/del", "V_DOMAIN_Functionality_comment","V_REGION_insertions", "V_REGION_deletions", 
                  "J_GENE_and_allele", "J_REGION_score", "J_REGION_identity_%", "J_GENE_and_allele_comment",  "JUNCTION_frame",
                  "Analysed_sequence_length","Locus")
  colnames(IMGTsum) <- c("sequence_number", "barcode", "contig", "sample", "functionality", "v_call", "v_score", "v_identity", "potential_indels", "v_comment", 
                         "v_ins", "v_del", "j_call", "j_score", "j_identity", "j_comment", "junction_frame", "analysed_length", "locus")
  
  IMGTseq <- imgtdfs[["2_IMGT-gapped-nt-sequences.txt"]] %>%
    dplyr::select("Sequence_number", "V_D_J_REGION", "V_J_REGION", "V_REGION") %>%
    dplyr::mutate(sequence_alignment = coalesce(as.character(V_D_J_REGION), as.character(V_J_REGION))) %>%
    dplyr::select(-"V_D_J_REGION", -"V_J_REGION")
  colnames(IMGTseq) <- c("sequence_number", "v_alignment", "sequence_alignment")
  
  IMGTnt <- imgtdfs[["3_Nt-sequences.txt"]] %>%
    dplyr::select("Sequence_number", "V_REGION_partial_5prime_missing_nt_nb_", "V_REGION_uncertain_nt_nb", "J_REGION_partial_3prime_missing_nt_nb",
                  "V_D_J_REGION", "V_J_REGION", "D_J_REGION", "V_REGION", "J_REGION", "JUNCTION") %>%
    dplyr::mutate(sequence = coalesce(as.character(V_D_J_REGION), as.character(V_J_REGION))) %>%
    dplyr::select(-"V_D_J_REGION", -"V_J_REGION")
  colnames(IMGTnt) <- c("sequence_number", "v_partial_missing", "v_uncertain_nt", "j_partial_missing", "d_j_region", "v_region","j_region","junction_nt", "sequence")
  
  IMGTaa <- imgtdfs[["4_IMGT-gapped-AA-sequences.txt"]] %>%
    dplyr::select("Sequence_number", "V_D_J_REGION", "V_J_REGION", "V_REGION") %>%
    dplyr::mutate(sequence = coalesce(as.character(V_D_J_REGION), as.character(V_J_REGION))) %>%
    dplyr::select(-"V_D_J_REGION", -"V_J_REGION")
  colnames(IMGTaa) <- c("sequence_number", "v_alignment_aa", "sequence_alignment_aa")
  
  IMGTAA <- imgtdfs[["5_AA-sequences.txt"]] %>%
    dplyr::select("Sequence_number", "V_D_J_REGION", "V_J_REGION", "JUNCTION") %>%
    dplyr::mutate(sequence = coalesce(as.character(V_D_J_REGION), as.character(V_J_REGION))) %>%
    dplyr::select(-"V_D_J_REGION", -"V_J_REGION") %>%
    dplyr::mutate(hydrophobicity = Peptides::hydrophobicity(imgtdfs[[6]]$`JUNCTION`, scale = "KyteDoolittle"))
  colnames(IMGTAA) <- c("sequence_number","junction_aa", "sequence_aa", "junction_hydrophobicity")
  
  IMGTjunc <- imgtdfs[["6_Junction.txt"]] %>%
    dplyr::select("Sequence_number", "pI") 
  colnames(IMGTjunc) <- c("sequence_number", "junction_pi")
  
  IMGTmut <- imgtdfs[["8_V-REGION-nt-mutation-statistics.txt"]] %>%
    dplyr::select("Sequence_number", "V_REGION_Nb_of_mutations", "V_REGION_Nb_of_silent_mutations", "V_REGION_Nb_of_nonsilent_mutations") %>%
    dplyr::mutate(V_REGION_Nb_of_mutations = as.numeric(gsub(" .*", "", V_REGION_Nb_of_mutations))) %>%
    dplyr::mutate(V_REGION_Nb_of_silent_mutations = as.numeric(gsub(" .*", "", V_REGION_Nb_of_silent_mutations))) %>%
    dplyr::mutate(V_REGION_Nb_of_nonsilent_mutations = as.numeric(gsub(" .*", "", V_REGION_Nb_of_nonsilent_mutations)))
  
  colnames(IMGTmut) <- c("sequence_number", "total_mut", "s_mut", "r_mut")
  IMGTmut
  
  #IMGTmutstat <- ParseIMGTmut(imgtdfs[["7_V-REGION-mutation-and-AA-change-table.txt"]])
  #colnames(IMGTmutstat) <- c("sequence_number", "fr1_mutlist", "cdr1_mutlist", "fr2_mutlist", "cdr2_mutlist", "fr3_mutlist")
  
  IMGTout <- IMGTsum %>%
    left_join(IMGTseq, by = "sequence_number") %>%
    left_join(IMGTnt, by = "sequence_number" ) %>%
    left_join(IMGTjunc, by = "sequence_number") %>%
    left_join(IMGTaa, by = "sequence_number") %>%
    left_join(IMGTAA, by = "sequence_number") %>%
    left_join(IMGTmut, by = "sequence_number") %>%
    #left_join(IMGTmutstat, by = "sequence_number") %>%
    filter(!grepl(" P", v_call))
  
  return(IMGTout)
}


#' @export
GetIMGT <- function(dir){
  IMGTdfs <- list()
  files <- list.files(path = dir, pattern="*.txt", full.names=TRUE, recursive=FALSE) 
  files <- files[!grepl("README.txt|11_Parameters.txt", files)]
  for (i in 1:length(files)){
    nam <- paste(gsub(".*/", "", files))
    IMGTdfs[[paste0(nam[i])]] <- data.table::as.data.table(readIMGT(files[i]))
    IMGTdfs[[i]] <- IMGTdfs[[i]] %>%
      tidyr::separate(Sequence_ID, c("barcode", "contig", "sample"), "-") %>%
      dplyr::mutate(Locus = ifelse(grepl("IGH", V_GENE_and_allele), "Heavy", "Light"))
  }
  IMGTout <- compileIMGTdfs(IMGTdfs)
  IMGTout <- IMGTout %>% mutate(vfam = gsub("\\-.*", "", v_call))
  
  #Parsing locus
  cat("Parsing locus")
  if (any(IMGTout$locus == "Light")) {
    Heavy <- data.table::as.data.table(IMGTout %>% filter(locus == "Heavy"))
    Light <- data.table::as.data.table(IMGTout %>% filter(locus == "Light"))
    locusout <- list(Heavy, Light)
    names(locusout) <- c("Heavy", "Light")
  } else {
    locusout <- IMGTout %>% filter(locus == "Heavy")
  }
  return(locusout)
}

