#' @export find_best_trim_lengths
#' @import dada2
#' @import ShortRead

find_best_trim_lengths <- function(
    directory=".",
    forward.pattern="_R1.fq.gz", # pattern for forward files, if single-end reverse reads are to be analyzed this should be a pattern for reverse reads
    reverse.pattern="_R2.fq.gz", # optional pattern for reverse files
    quality.threshold=27, # minimal acceptable average base quality, 27 means ~ 1 error in 500 bases
    use.mean=F, # reads are trimmed if mean per-cycle score falls below this value
    max.percent.dropped=10, # reads are trimmed at such a point that max. x% is dropped, i.e. min 100 - x% trimmed reads must have all scores over quality.threshold
    amplicon.length=450,
    min.overlap=50,
    min.length=250, # mininmum length for single-end reads
    nfiles=10
){

  # function to get fraction of reads with score better than quality.threshold
  get_quant <- function(xx, q) {
    (cumsum(xx)/sum(xx))[[q]]
  }
  out <- c(0, 0) # no trimming

    fl <- list.files(directory, forward.pattern)
    if(length(fl) > nfiles){
      fl <- fl[1:nfiles]
    }
    qaF <- ShortRead::qa(fl, n=10000)
    qaFdf <- as.data.frame(qaF[["perCycle"]]$quality)
    meansF <- rowsum(qaFdf$Score * qaFdf$Count, qaFdf$Cycle)/rowsum(qaFdf$Count, qaFdf$Cycle)
    qaFdf_aggregated <- aggregate(Count ~ Score * Cycle, data=qaFdf, FUN='sum')
    qF <- by(qaFdf_aggregated, qaFdf_aggregated$Cycle, function(x) sum(x[x$Score < quality.threshold, ]$Count)/sum(x$Count))
    cumF <- by(qaFdf, qaFdf$Cycle, function(foo) sum(foo$Count), simplify = TRUE)
    statdf <- as.data.frame(cbind(qF, meansF, cumF))
    colnames(statdf) <- c("qF","meansF", "cumF")
    statdf$cycle <-as.numeric(rownames(statdf))
    statdfF <- "statdfF"
    assign(statdfF, statdf, envir=parent.frame())
    if(use.mean){
      statdfcut <- statdf[statdf$meansF < quality.threshold, ]
      cutoffF <- min(statdfcut$cycle)
    }else{
      statdfcut <- statdf[statdf$qF > max.percent.dropped/100, ]
      cutoffF <- min(statdfcut$cycle)
      message
    }

    if(!is.null(reverse.pattern)){
      fr <- list.files(directory, reverse.pattern)
      if(length(fr) > nfiles){
        fr <- fr[1:nfiles]
      }
      qaR <- ShortRead::qa(fr, n=10000)
      qaRdf <- as.data.frame(qaR[["perCycle"]]$quality)
      meansR <- rowsum(qaRdf$Score * qaRdf$Count, qaRdf$Cycle)/rowsum(qaRdf$Count, qaRdf$Cycle)
      qaRdf_aggregated <- aggregate(Count ~ Score * Cycle, data=qaRdf, FUN='sum')
      qR <- by(qaRdf_aggregated, qaRdf_aggregated$Cycle, function(x) sum(x[x$Score < quality.threshold, ]$Count)/sum(x$Count))
      cumR <- by(qaRdf, qaRdf$Cycle, function(foo) sum(foo$Count), simplify = TRUE)
      statdf <- as.data.frame(cbind(qR, meansR, cumR))
      colnames(statdf) <- c("qR","meansR", "cumR")
      statdf$cycle <-as.numeric(rownames(statdf))
      statdfR <- "statdfR"
      assign(statdfR, statdf, envir=parent.frame())
      if(use.mean){
        statdfcut <- statdf[statdf$meansR < quality.threshold, ]
        cutoffR <- min(statdfcut$cycle)
      }else{
        statdfcut <- statdf[statdf$qR > max.percent.dropped/100, ]
        cutoffR <- min(statdfcut$cycle)
      }
      if(amplicon.length > cutoffF + cutoffR - min.overlap){
        stop("find_best_trim_lengths: trimmed lengths do not ensure enough overlap to merge paired reads")
      }else{
        out <- c(cutoffF, cutoffR)
      }
    }else{ # we're processing single-end data
      if(cutoffF < min.lenght){
        stop("find_best_trim_lengths: acceptable trimmed length was lower than min.length threshold of ", min.length)
      }else{
        out <- c(cutoffF)
      }
    }

  return(out)
}
