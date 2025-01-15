#' @export fasta.file
#'

fasta.file <- function(
    experiment=NULL,
    feature=NULL,
    d=NULL,
    rarefied=NULL){
  out <- NULL
  out <- mgmisc1::extract(experiment=experiment, feature=feature, d=d, rarefied=rarefied, what='fasta.file')
  return(out)
}
