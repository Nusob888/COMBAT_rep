##Functions for generating matrix of BCR metrics
#' @import mosaic
#' @import DescTools
#' @import parallel
#' @import data.table
#' @import magrittr

#Calculate shannons ----
#Shannon function - this is required for the actual shannon calculation 
shannon <- function(data) {
  counts <- table(data)
  p.i <- counts/sum(counts)
  H <- (-1) * sum(p.i * log(p.i))
  return(H)
}
#Function to counter Gini behaviour when resampling returns a single clone.. if a vector length =1, Gini generates NaN instead of 1 (100% inequality)----
custom_resample <- function(data, subsample=subsample){
  counts<- table(mosaic::resample(data, subsample, replace=FALSE))
  if(length(counts)==1){
    custom_count <- c(counts,0)
    return(custom_count)
  }else{
    return(counts)
  }
}
##Bootstrap shannon function
boot.shannon <- function(input, repertoire_id, cell_subset, clone_id, B = 1000, subsample= subsample,conf.level = 0.95){
  if (is.null(cell_subset)){
  data <- input[[clone_id]]
  ID <- input[[repertoire_id]] %>% unique #this is just to set the IDs later. If you add a new variable to your group split above, you will need to add the same variable at the end of here e.g. if I add cluster_id then that needs to go after "Source"
  options(`mosaic:parallelMessage` = FALSE)
  shannon.boot <- mosaic::do(B, parallel=TRUE) * shannon(mosaic::resample(data, subsample, replace=FALSE)) #Here we subsample with replacement and calculate a shannons for those 30 subsamples. This is repeated x number of B
  names(shannon.boot) <- "result"
  boot.mean <- mosaic::mean(~result, data = shannon.boot)
  boot.se <- mosaic::sd(~result, data = shannon.boot)
  boot.median <- mosaic::median(~result, data = shannon.boot)
  g <- as.numeric(1 - conf.level)
  boot.lci <- mosaic::qdata(~result, g/2,data = shannon.boot)
  boot.uci <- mosaic::qdata(~result,1-g, data = shannon.boot)
  boot.stats <- data.table::data.table(mean.shannon=boot.mean, se=boot.se, median.shannon=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
  boot.stats <- cbind(ID, boot.stats)
  return(boot.stats)
  } else {
    data <- input[[clone_id]]
    ID <- input[,c(repertoire_id, cell_subset)] %>% unique #this is just to set the IDs later. If you add a new variable to your group split above, you will need to add the same variable at the end of here e.g. if I add cluster_id then that needs to go after "Source"
    options(`mosaic:parallelMessage` = FALSE)
    shannon.boot <- mosaic::do(B) * shannon(mosaic::resample(data, subsample, replace=FALSE)) #Here we subsample with replacement and calculate a shannons for those 30 subsamples. This is repeated x number of B
    names(shannon.boot) <- "result"
    boot.mean <- mosaic::mean(~result, data = shannon.boot)
    boot.se <- mosaic::sd(~result, data = shannon.boot)
    boot.median <- mosaic::median(~result, data = shannon.boot)
    g <- as.numeric(1 - conf.level)
    boot.lci <- mosaic::qdata(~result, g/2,data = shannon.boot)
    boot.uci <- mosaic::qdata(~result,1-g, data = shannon.boot)
    boot.stats <- data.table::data.table(mean.shannon=boot.mean, se=boot.se, median.shannon=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
    boot.stats <- cbind(ID, boot.stats)
    return(boot.stats)
  }
}

##CEI----
getCEI <- function(input, repertoire_id, cell_subset, clone_id, B = 1000, subsample = subsample, conf.level = 0.95) {
  if (is.null(cell_subset)){
    data <- input[[clone_id]]
    ID <- input[[repertoire_id]] %>% unique #this is just to set the IDs later. If you add a new variable to your group split above, you will need to add the same variable at the end of here e.g. if I add cluster_id then that needs to go after "Source"
    options(`mosaic:parallelMessage` = FALSE)
  gini.boot <- mosaic::do(B) * DescTools::Gini(custom_resample(data, subsample=subsample))
  names(gini.boot) <- "result"
  #As Gini can only accept a vector of more than two, when we only detect one clone, it creates a NaN by coercian. 
  #To fix this, we will manually reset NaNs to 1
  boot.mean <- mosaic::mean(~result, data = gini.boot)
  boot.se <- mosaic::sd(~result, data = gini.boot)
  boot.median <- mosaic::median(~result, data = gini.boot)
  g <- as.numeric(1 - conf.level)
  boot.lci <- mosaic::qdata(~result, g/2,data = gini.boot)
  boot.uci <- mosaic::qdata(~result,1-g, data = gini.boot)
  boot.stats <- data.table::data.table(mean.CEI=boot.mean, se=boot.se, median.CEI=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
  boot.stats <- cbind(ID, boot.stats)
  return(boot.stats)
  } else {
    data <- input[[clone_id]]
    ID <- input[,c(repertoire_id, cell_subset)] %>% unique #this is just to set the IDs later. If you add a new variable to your group split above, you will need to add the same variable at the end of here e.g. if I add cluster_id then that needs to go after "Source"
    options(`mosaic:parallelMessage` = FALSE)
    gini.boot <- mosaic::do(B) * DescTools::Gini(custom_resample(data, subsample=subsample))
    names(gini.boot) <- "result"
    #As Gini can only accept a vector of more than two, when we only detect one clone, it creates a NaN by coercian. 
    #To fix this, we will manually reset NaNs to 1
    boot.mean <- mosaic::mean(~result, data = gini.boot)
    boot.se <- mosaic::sd(~result, data = gini.boot)
    boot.median <- mosaic::median(~result, data = gini.boot)
    g <- as.numeric(1 - conf.level)
    boot.lci <- mosaic::qdata(~result, g/2,data = gini.boot)
    boot.uci <- mosaic::qdata(~result,1-g, data = gini.boot)
    boot.stats <- data.table::data.table(mean.CEI=boot.mean, se=boot.se, median.CEI=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
    boot.stats <- cbind(ID, boot.stats)
    return(boot.stats)
  }
}

##CDI

cdi_resample <- function(data, subsample=subsample, clone_id=clone_id){
  var1 <-  enquo(clone_id)
  sub <- mosaic::resample(data, subsample, replace=FALSE)
  table <- table(unique(sub)[[clone_id]])[as.character(sub[[clone_id]])]
  if(length(table)==1){
    custom_count <- c(table,0)
    return(custom_count)
  }else{
  return(table)
  }
}

getCDI <- function(input, repertoire_id, cell_subset, clone_id, B = 1000, subsample = subsample, conf.level = 0.95) {
  if (is.null(cell_subset)){
    data <- input %>% select(all_of(clone_id), sequence) 
    ID <- input[[repertoire_id]] %>% unique #this is just to set the IDs later. If you add a new variable to your group split above, you will need to add the same variable at the end of here e.g. if I add cluster_id then that needs to go after "Source"
    options(`mosaic:parallelMessage` = FALSE)
    gini.boot <- mosaic::do(B) * DescTools::Gini(cdi_resample(data=data, clone_id=clone_id, subsample=subsample))
    names(gini.boot) <- "result"
    #As Gini can only accept a vector of more than two, when we only detect one clone, it creates a NaN by coercian. 
    #To fix this, we will manually reset NaNs to 1
    boot.mean <- mosaic::mean(~result, data = gini.boot)
    boot.se <- mosaic::sd(~result, data = gini.boot)
    boot.median <- mosaic::median(~result, data = gini.boot)
    g <- as.numeric(1 - conf.level)
    boot.lci <- mosaic::qdata(~result, g/2,data = gini.boot)
    boot.uci <- mosaic::qdata(~result,1-g, data = gini.boot)
    boot.stats <- data.table::data.table(mean.CDI=boot.mean, se=boot.se, median.CDI=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
    boot.stats <- cbind(ID, boot.stats)
    return(boot.stats)
  } else {
    data <- input %>% select(all_of(clone_id), sequence) 
    ID <- input[,c(repertoire_id, cell_subset)] %>% unique #this is just to set the IDs later. If you add a new variable to your group split above, you will need to add the same variable at the end of here e.g. if I add cluster_id then that needs to go after "Source"
    options(`mosaic:parallelMessage` = FALSE)
    gini.boot <- mosaic::do(B) * DescTools::Gini(cdi_resample(data=data, clone_id=clone_id, subsample=subsample))
    names(gini.boot) <- "result"
    #As Gini can only accept a vector of more than two, when we only detect one clone, it creates a NaN by coercian. 
    #To fix this, we will manually reset NaNs to 1
    boot.mean <- mosaic::mean(~result, data = gini.boot)
    boot.se <- mosaic::sd(~result, data = gini.boot)
    boot.median <- mosaic::median(~result, data = gini.boot)
    g <- as.numeric(1 - conf.level)
    boot.lci <- mosaic::qdata(~result, g/2,data = gini.boot)
    boot.uci <- mosaic::qdata(~result,1-g, data = gini.boot)
    boot.stats <- data.table::data.table(mean.CDI=boot.mean, se=boot.se, median.CDI=boot.median, lower_ci=boot.lci, upper_ci=boot.uci)
    boot.stats <- cbind(ID, boot.stats)
    return(boot.stats)
  }
}




#Bootstrapping function - wrapper to bootstrap ----
# input = dataframe
# B = number of iterations, preset at 1000 
# conf.level = preset as a standard 0.95 CI
# subsample = subsampling size
#' @export
get_bcr_metrics <- function(input, method = NULL, repertoire_id=NULL, cell_subset=NULL, clone_id=NULL, B = 1000, subsample= 17,conf.level = 0.95) {
  numCores <- parallel::detectCores()-1
  if (inherits(input, "list") == FALSE){
    warning("input is not a list of data frames, please input list of dataframes groupsplit by repertoire_id")
  }else if (is.null(clone_id)) {
    print("Please input clone_id")
  } else if (is.null(repertoire_id)) {
      print("Please state repertoire_id")
    } else if (is.null(method)) {
       print("Please input method of choice from 'shannon', 'CDI', 'CEI' or 'all' to return a list of dataframes") 
    } else {
      if (method=="shannon"){
          out <- parallel::mclapply(input,function(x){boot.shannon(input=x, repertoire_id=repertoire_id, cell_subset=cell_subset, clone_id=clone_id,subsample=subsample)}, mc.cores=numCores)
          out <- data.table::rbindlist(out) 
          return(out)
      }
      
      if (method=="CEI"){
        out <- parallel::mclapply(input,function(x){getCEI(input=x, repertoire_id=repertoire_id, cell_subset=cell_subset, clone_id=clone_id,subsample=subsample)}, mc.cores=numCores)
        out <- data.table::rbindlist(out) 
      return(out)
      }
      
      if (method=="CDI"){
        out <- parallel::mclapply(input,function(x){getCDI(input=x, repertoire_id=repertoire_id, cell_subset=cell_subset, clone_id=clone_id,subsample=subsample)}, mc.cores=numCores)
        out <- data.table::rbindlist(out) 
        return(out)
      }
      
      if (method=="all"){
        shannon <- parallel::mclapply(input,function(x){boot.shannon(input=x, repertoire_id=repertoire_id, cell_subset=cell_subset, clone_id=clone_id,subsample=subsample)}, mc.cores=numCores)
        shannon <- data.table::rbindlist(shannon) 
        
        CEI <- parallel::mclapply(input,function(x){getCEI(input=x, repertoire_id=repertoire_id, cell_subset=cell_subset, clone_id=clone_id,subsample=subsample)}, mc.cores=numCores)
        CEI <- data.table::rbindlist(CEI) 
        
        CDI <- parallel::mclapply(input,function(x){getCDI(input=x, repertoire_id=repertoire_id, cell_subset=cell_subset, clone_id=clone_id,subsample=subsample)}, mc.cores=numCores)
        CDI <- data.table::rbindlist(CDI) 
        
        out <- list(shannon, CEI, CDI)
        names(out) <- c("shannon", "CEI", "CDI")
        return(out)
      }
    }
  
}







