#' @export generateOneSlice
#'

generateOneSlice <- function(
    experiment=NULL,
    variable=NULL,
    l=NULL){
  sdata <- mgmisc1::sampledata(experiment=experiment, feature="original.data")
  split_sdata <- sdata[ sdata[[variable]] == l, ]
  if(nrow(split_sdata) < 3){
    outlist <- list()
    message("generateOneSlice: slice ", l, " has less than 3 rows, skipping it")
    return(outlist)
  }
  outlist <- mgmisc1::new_mg()
  outlist$metadata <- experiment$metadata
  for( feature in c('original.data', 'ASVs', 'OTUs', 'KOs', 'pathways')){
    message("generateOneSlice: slicing ", feature)
    if(feature != 'OTUs'){
      f <- mgmisc1::extractFeatureSet(experiment=experiment, feature=feature)
      outlist[[feature]] <- mgmisc1::sliceFeatureSet(set=f, split_sampledata=split_sdata, split_featureannot=f$featureannot)
      if(feature %in% c('KOs', 'pathways')){
        outlist[[feature]]$descriptions <- outlist[[feature]]$taxonomy
        outlist[[feature]]$taxonomy <- NULL
      }
      if(!is.null(experiment[[feature]]$rarefied)){
        f <- mgmisc1::extractFeatureSet(experiment=experiment, feature=feature, rarefied=T)
        outlist[[feature]]$rarefied <- mgmisc1::sliceFeatureSet(set=f, split_sampledata=split_sdata, split_featureannot=f$featureannot)
        if(feature %in% c('KOs', 'pathways')){
          outlist[[feature]]$rarefied$descriptions <- outlist[[feature]]$rarefied$taxonomy
          outlist[[feature]]$rarefied$taxonomy <- NULL
        }
      }
    }else{
      for( d in names(experiment$OTUs$dists)){
        message("generateOneSlice: slicing ", d, " OTUs")
        f <- extractFeatureSet(experiment=experiment, feature=feature, d=d)
        outlist[[feature]]$dists[[d]] <- mgmisc1::sliceFeatureSet(set=f, split_sampledata=split_sdata, split_featureannot=f$featureannot)
        if(!is.null(experiment[[feature]]$dists[[d]]$rarefied)){
          f <- mgmisc1::extractFeatureSet(experiment=experiment, feature=feature, d=d, rarefied=T)
          outlist[[feature]]$dists[[d]]$rarefied <- mgmisc1::sliceFeatureSet(set=f, split_sampledata=split_sdata, split_featureannot=f$featureannot)
        }
      }
    }
  }

  return(outlist)
}
