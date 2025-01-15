#' @export extract
#'


extract <- function(
    experiment=NULL,
    feature=NULL,
    diss=0.03, # dissimilarity threshold, needed for OTUs
    rarefied=F,
    what=NULL
){
  if(!(feature %in% c('ASVs', 'OTUs', 'KOs', 'pathways', 'original.data'))){
    stop("extract: feature to extract from must be on of 'ASVs', 'OTUs', 'KOs', 'pathways', 'original.data'")
  }
  if(!(what %in% c('featuretab', 'sampledata', 'tree', 'tree.file', 'featureannot', 'descriptions', 'fasta', 'fasta.file', 'count_table', 'count.file'))){
    stop("extract: what to extract must be one of 'featuretab', 'sampledata', 'tree', 'tree.file', 'featureannot', 'descriptions', 'fasta', 'fasta.file', 'count_table', 'count.file'")
  }
  if(is.null(experiment[[feature]])){
    message("extract: nothing to extract from")
    return(NULL)
  }

  diss <- as.character(diss)
  out <- NULL
  if(feature == 'OTUs' & !is.null(diss)){
    if(!rarefied){
      out <- experiment[[feature]]$dists[[diss]][[what]]
    }else{
      out <- experiment[[feature]]$dists[[diss]]$rarefied.data[[what]]
    }
  }else if( feature != 'OTUs'){
    if(!rarefied){
      out <- experiment[[feature]][[what]]
    }else{
      out <- experiment[[feature]]$rarefied.data[[what]]
    }
  }
  else{
    stop("extract: if feature to be extracted is 'OTUs', dissimilarity level (diss) must be also given")
  }
    if(feature == "original.data" & what == "featuretab"){
      message("extract: extracting seqtab")
      out <- experiment$original.data$seqtab
    }

  return(out)
}
