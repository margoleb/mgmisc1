#' @export featureNames

featureNames <- function(
  experiment = NULL,
  feature = NULL,
  rarefied = FALSE,
  diss = c(0.03)
  ){
  if(feature == "OTUs"){
    fnames <- list()
  }else{
    fnames <- c()
  }

  if(!is.null(experiment$feature)){
    if(feature != "OTUs"){
      if(rarefied){
        fnames <- colnames(experiment[[feature]]$rarefied$featuretab)
      }else{
        fnames <- colnames(experirment[[feature]]$featuretab)
      }
    }else{
      for(d in diss){
        if(rarefied){
          fnames[[d]] <- colnames(experiment$OTUs$rarefied$dists[[d]]$featuretab)
        }else{
          fnames[[d]] <- colnames(experiment$OTUs$dists[[d]]$featuretab)
        }
      }
    }
  }else{
    message("featureNames: no features, empty vector returned")
  }

  return(fnames)
}

