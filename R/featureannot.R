#' @export featureannot
#'

featureannot <- function(
    experiment=NULL,
    feature=NULL,
    d=NULL,
    rarefied=F){
  out <- NULL
  out <- mgmisc1::extract(experiment=experiment, feature=feature, d=d, rarefied=rarefied, what='featureannot')
  return(out)
}
