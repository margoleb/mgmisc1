#' @export sampleNames


sampleNames <- function(
    experiment = NULL,
    feature = NULL,
    rarefied = FALSE,
    diss = c("0.03")
  ){
  if(feature == "OTUs"){
    snames <- list()
  }else{
    snames <- c()
  }

  if(!is.null(experiment[[feature]])){
    if(feature != "OTUs"){
      if(rarefied){
        snames <- rownames(experiment[[feature]]$rarefied.data$featuretab)
      }else{
        snames <- rownames(experiment[[feature]]$featuretab)
      }
    }else{
      for(d in diss){
        if(rarefied){
          snames[[d]] <- rownames(experiment$OTUs$rarefied.data$dists[[d]]$featuretab)
        }else{
          snames[[d]] <- rownames(experiment$OTUs$dists[[d]]$featuretab)
        }
      }
    }
  }else{
    message("sampleNames: no samples, empty vector returned")
  }

  return(snames)

}

