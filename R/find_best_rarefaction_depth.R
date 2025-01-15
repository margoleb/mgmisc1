#' @export find_best_rarefaction_depth


find_best_rarefaction_depth <- function(
  featuretab=NULL,
  sample.weight=3000, # how many times a sample is more important than a sequence/read
  min.depth=500,
  step=10, # interval at which calculations are performed
  max.samples.discarded=10 # maximal percentage of the number of samples discarded
  ){
    max.score = 0
    best.depth = min.depth
    best.nsamples = 1
    read_numbers <- rowSums(featuretab)
    for(depth in seq(min.depth, max(read_numbers), step)){
      nsamples <- sum(read_numbers > depth)
      score <- nsamples * (sample.weight + depth)
      if(score >= max.score){
        max.score = score
        best.depth = depth
        best.nsamples = nsamples
      }
    }
    nsamples.discarded <- length(read_numbers) - best.nsamples
    percentage.discarded <- 100 * nsamples.discarded/length(read_numbers)
    if(percentage.discarded > max.samples.discarded ){
      message("Finding rarefaction depth satisfying max percentage of samples discarded was impossible, try decreasing sample.weight and/or min.depth and increase max.samples.discarded ")
      best.depth = 0
    }

  return(best.depth)
}
