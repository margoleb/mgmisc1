#' @export sliceFeatureSet
#'

sliceFeatureSet <- function(
    set = NULL, # list comprising 'featuretab', 'sampledata', 'featureannot', 'tree', 'fasta', 'fasta.file', 'count_table', and 'count.file' components, obtained with extractFeatureSet
    split_sampledata = NULL,
    split_featureannot = NULL
){
  outlist <- list()
  if(is.null(split_featureannot)){
    split_featureannot <- set$featureannot
  }
  outlist <- mgmisc1::makeConsistent(featuretab=set$featuretab,
                                     sampledata=split_sampledata,
                                     featureannot=split_featureannot,
                                     tree=set[["tree"]],
                                     tree.file=set[["tree.file"]],
                                     fasta=set[["fasta"]],
                                     fasta.file=set[["fasta.file"]],
                                     count_table=set$count_table,
                                     count.file=set$count.file)
  return(outlist)
}
