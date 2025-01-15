#' @export extractFeatureSet
#'

extractFeatureSet <- function(
    experiment=NULL,
    feature=NULL,
    rarefied=F,
    d=NULL # distance threshold needed if extracting OTUs
){
  outlist <- list(featuretab=NULL, featureannot=NULL, sampledata=NULL, fasta=NULL, fasta.file="", count_table=NULL, count.file="", tree=NULL, tree.file="")
  if(class(experiment) != "mg"){
    stop("extractFeatureSet: an experiment object must be of class 'mg'")
  }
  # If there is no requested feature in the experiment object, return NULL, it is needed when analyzing experiments to indicate a need to generate some data.
  if(is.null(experiment[[feature]])){
    message("extractFeatureSet: no feature named ", feature, " in analyzed object")
    return(NULL)
  }

  if(feature == "OTUs" && length(experiment$OTUs$dists[[d]]) == 0){
    message("extractFeatureSet: no OTUs at ", d, " dissimilarity in analyzed object")
    return(NULL)
  }

  if(feature %in% c("ASVs", "OTUs")){
    cat("extractFeatureSet: extracting fasta... ")
    outlist$fasta <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="fasta")
    message("done")
    cat("extractFeatureSet: extracting name of a fasta file... ")
    outlist$fasta.file <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="fasta.file")
    message(outlist$fasta.file, ", done")
    if(!(is.null(outlist$fasta))){
      if(!file.exists(outlist$fasta.file)){
        if(outlist$fasta.file != ""){
          write.table(outlist$fasta, outlist$fasta.file, quote=F, sep="", col.names=F)
        }
        else{
          outlist$fasta.file <- "temp.fasta"
          write.table(outlist$fasta, outlist$fasta.file, quote=F, sep="", col.names=F)
        }
      }
    }
    cat("extractFeatureSet: extracting count table... ")
    outlist$count_table <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="count_table")
    message("done")
    cat("extractFeatureSet: extracting count file... ")
    outlist$count.file <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="count.file")
    message(outlist$count.file, ", done")
    if(!is.null(outlist$count_table) && !file.exists(outlist$count.file)){
      write.table(outlist$count_table, outlist$count.file, quote=F, sep="\t", row.names=F)
    }
    cat("extractFeatureSet: extracting tree... ")
    outlist$tree <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="tree")
    message("done")
    cat("extractFeatureSet: extracting tree file... ")
    outlist$tree.file <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="tree.file")
    message(outlist$tree.file, ", done")
    if(!is.null(outlist$tree) && !file.exists(outlist$tree.file)){
      ape::write.tree(outlist$tree, outlist$tree.file)
    }
  }

  message("extractFeatureSet: extracting feature table")
  outlist$featuretab <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="featuretab")
  message("extractFeatureSet: extracting feature annotations")
  outlist$featureannot <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="featureannot")
  message("extractFeatureSet: extracting sample data")
  outlist$sampledata <- mgmisc1::extract(experiment=experiment, feature=feature, rarefied=rarefied, diss=d, what="sampledata")

  if(is.null(outlist$featuretab) || is.null(outlist$featureannot) || is.null(outlist$sampledata)){
    message("extractFeatureSet: one of featuretab, featureannot or sampledata was null")
    return(NULL)
  }

  return(outlist)
}


