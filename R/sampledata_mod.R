#' @export "sampledata<-.mg"
#'
# Modifies 'sampledata' slot of all features but 'original.data'

`sampledata<-.mg` <- function(
    experiment,
    ...,
    sdata){
  sd <- sampledata(experiment=experiment, feature="ASVs", rarefied=F)
  # are samplenames and their order the same in the original and new table?
  if(order(rownames(sd)) != order(rownames(sdata))){
    stop("sampledata: sets of samples in the original sampledata and the new one differ")
  }
  # non-rarefied data
  sdata <- sdata[ order(rownames(sd)), ]
  for(f in c('ASVs', 'OTUs', 'KOs', 'pathways')){
    if(f == 'OTUs'){
      for(d in names(experiment$OTUs$dists)){
        experiment$OTUs$dists[[d]]$sampledata <- sdata
      }
    }else{
      experiment[[f]]$sampledata <- sdata
    }
  }

  # rarefied data
  sd <- sampledata(experiment, feature="ASVs", rarefied=T)
  if(!is.null(sd)){
     sdata <- sdata[ rownames(sdata) %in% rownames(sd), ]
     sdata <- sdata[ rownames(sd), ]
     for(f in c('ASVs', 'OTUs', 'KOs', 'pathways')){
       if(f == 'OTUs'){
         for(d in names(experiment$OTUs$dists)){
           experiment$OTUs$dists[[d]]$rarefied.data$sampledata <- sdata
         }
       }else{
         experiment[[f]]$rarefied.data$sampledata <- sdata
       }
     }
  }
}
