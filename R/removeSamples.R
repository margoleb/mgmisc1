#' @export removeSamples

`%nin%` = Negate(`%in%`)

removeSamples <- function(
    experiment = NULL,
    samples = c())
{
  if(is.null(experiment)){
    stop("removeSamples: no experiment to work on")
  }
  samples <- paste0(samples, collapse='", "')
  e <- paste0('samplename %nin% c("', samples, '")')
  result <- subsetExperiment(experiment=experiment, what="samples", condition=e )

  return(result)
}
