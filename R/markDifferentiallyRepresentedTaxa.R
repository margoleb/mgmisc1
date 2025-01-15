#' @export markDifferentiallyRepresentedTaxa
#' @import ggplot2 ggtext
#'


markDifferentiallyRepresentedTaxa <- function(experiment=NULL,
                                              featurename="ASVs" # must be either "ASVs" or "OTUs"
                                              ){

  if(is.null(experiment[[featurename]]$analyzes$taxonomy$taxtabs) | !is.null(esperiment[[featurename]]$analyzes$diff.features)){

  }
  p <- p + scale_fill_discrete(labels=c()) # labels with mark-ups
  return(experiment)
}
