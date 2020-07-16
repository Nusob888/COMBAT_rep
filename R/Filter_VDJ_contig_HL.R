filterHeavycontigs<- function(test){
  empty <- data.frame(barcode_plex=NA,
                      IMGTchain=NA, 
                      umis=NA, 
                      umi_ranks=NA, 
                      VDJ_group=NA, 
                      c_gene=NA)
  
  #Parse chains
  Heavy <- test %>% filter(grepl("IGH", IMGTchain))
  Light <- test %>% filter(grepl("IGK|IGL", IMGTchain))

#Process Heavy chain
 if(exists("Heavy")){
  if(!any(Heavy$umi_ranks== "2")){
    x <- Heavy[c(Heavy$umi_ranks == "1"),] #select rank1 if rank2 ratio less than < 0.5 (rank1 umis is double of rank2)
  }else{ 
  for(i in 1:nrow(Heavy)){
     if (Heavy$umi_ranks[i] == "1"){#----
        if(Heavy$c_gene[i] == "IGHM"){#conditional selection of rank1 IGHMs
          if(Heavy[c(Heavy$umi_ranks== "1"), c("VDJ_group")] == Heavy[c(Heavy$umi_ranks== "2"), c("VDJ_group")]){
            x <- Heavy %>% 
              filter(grepl(paste(Heavy$VDJ_group[i]), VDJ_group)) %>% #filter on VJ and cdr3 matches
              filter(grepl("IGHM|IGHD", c_gene)) #filter VJ/cdr3 matches that are IGHM or IGHD, thereby for top IGHMs, keeping paired IGHD expression
          }else if((Heavy[c(Heavy$umi_ranks == "2"),c("umis")]/Heavy[c(Heavy$umi_ranks == "1"), c("umis")]) <0.4){
            x <- Heavy %>% 
              filter(grepl(paste(Heavy$VDJ_group[i]), VDJ_group)) %>% #filter on VJ and cdr3 matches
              filter(grepl("IGHM|IGHD", c_gene)) #filter VJ/cdr3 matches that are IGHM or IGHD, thereby for top IGHMs, keeping paired IGHD expression
          }else{
            x <- empty
            print(paste("No IGHM assignments possible for contig:", Heavy$barcode_plex[i], ", potential doublet", sep="")) 
          }
        }else if(test$c_gene[i] != "IGHM"){# conditional selection of non_IGHM rank1s
          if(Heavy[c(Heavy$umi_ranks== "1"), c("VDJ_group")] == Heavy[c(Heavy$umi_ranks== "2"), c("VDJ_group")]){
            x <- Heavy %>% 
              filter(grepl(paste(Heavy$VDJ_group[i]), VDJ_group))
          }else if((Heavy[c(Heavy$umi_ranks == "2"),c("umis")]/Heavy[c(Heavy$umi_ranks == "1"), c("umis")]) <0.4){
            x <- Heavy %>% 
              filter(grepl(paste(Heavy$VDJ_group[i]), VDJ_group)) 
          }else{
            x <- empty
            print(paste("No IGHG assignments possible for contig:", Heavy$barcode_plex[i], ", potential doublet", sep=""))  
          }
        }
     }
    }
  }
}

#Process light chain
if(exists("Light")){
  if(!any(Light$umi_ranks== "2")){
    y <- Light[c(Light$umi_ranks == "1"),] #select rank1 if rank2 ratio less than < 0.5 (rank1 umis is double of rank2)
  }else{
  for(i in 1:nrow(Light)){
    if (Light$umi_ranks[i] == "1"){
    if((Light[c(Light$umi_ranks == "2"),c("umis")]/Light[c(Light$umi_ranks == "1"), c("umis")]) <0.4){
      y <- Light %>% 
        filter(grepl(paste(Light$VDJ_group[i]), VDJ_group)) #filter VJ/cdr3 matches that are IGHM or IGHD, thereby for top IGHMs, keeping paired IGHD expression
    }else{
      y <- empty
      print(paste("No IGK/L assignments possible for contig:", Light$barcode_plex[i], ", potential doublet", sep=""))  
    }
    }
  }
    }
}

#Sort objects
if(exists("x") & exists("y")){
  xy <- mapply(c, x, y, SIMPLIFY=FALSE)
  return(xy)
}else if(exists("x") & isFALSE(exists("y"))){
  return(x)
}else if(exists("y") & isFALSE(exists("x")))
  return(y)
}