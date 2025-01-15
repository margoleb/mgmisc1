#' @export prepareBIOM


prepareBIOM <- function(
  featuretab=NULL,
  write = TRUE,
  filename = "bacteria.biom"
){
  result <- biomformat::make_biom(data = t(featuretab))
  if(write){
    biomformat::write_biom(result, filename)
  }

  return(result)
}
