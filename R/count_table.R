#' @export count_table
#'

count_table <- function(
    experiment=NULL,
    feature=NULL,
    d=NULL,
    rarefied=NULL){
  out <- NULL
  out <- mgmisc1::extract(experiment=experiment, feature=feature, d=d, rarefied=rarefied, what='count_table')
  return(out)
}
