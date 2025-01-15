#' @export tree
#'

tree <- function(
    experiment=NULL,
    feature=NULL,
    d=NULL,
    rarefied=NULL){
  out <- NULL
  out <- mgmisc1::extract(experiment=experiment, feature=feature, d=d, rarefied=rarefied, what='tree')
  return(out)
}
