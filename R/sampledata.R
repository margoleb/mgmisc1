#' @export sampledata
#'

sampledata <- function(
    experiment=NULL,
    feature=NULL,
    d=NULL,
    rarefied=F){
  out <- NULL
  out <- as.data.frame(mgmisc1::extract(experiment=experiment, feature=feature, d=d, rarefied=rarefied, what='sampledata'))
  return(out)
}

