#' denoiseFastq
#'
#' Filter, trim, denoise, and merge reads, identify and remove chimeras and construct sequence table.
#'
#' @param directory A path to a directory where fastq files are stored. Defaults to '.'.
#' @param novaseqLoess A function to use while denoising NovaSeq data if dada2 version is below 1.32. Current possibilities: loessErrfun_mod1, loessErrfun_mod2, loessErrfun_mod3, loessErrfun_mod4. Defaults to NULL (standard version works well with NovaSeq data as of dada2 v. 1.32).
#' @param save A logical indicating if resulting sequence table should be saved as an RDS file. Defaults to TRUE.
#' @param basename A string with which names of all generated files will begin. Defaults to 'bacteria'.
#' @param pattern A pattern common for all fastq files to be analyzed. Defaults to '_R1_001.fastq.gz'.
#' @param filtered A logical indicating if precomputed filtering results should be used. If TRUE there needs to be a 'filtered' directory in the dir where fastq files are stored and it must contain filtered fastq files. Defaults to FALSE.
#' @param error.rates A path to an RDS file containing a list consisting of errF and errR error matrices from dada2's learnErrors(). If not NULL precomputed error rates will be used.  Defaults to NULL.
#' @param sample.name.sep A string separating sample name from the rest of filename. Defaults to '_'.
#' @param truncLen A vector of two positive integers indicating where R1 and R2 should be truncated or 'auto'. If 'auto' optimal trimming positions are calculated (see find_best_trim_lengths()). For NovaSeq usually there is no need to truncate reads, so c(0,0) or c(250,250) should be used. Defaults to 'auto'.
#' @param pool A logical indicating if samples should be pooled which improves rare variants detection or 'pseudo'. Increases both memory consumption and CPU time, use with huge data is discouraged. Defaults to FALSE.
#' @param quality.threshold A positive integer indicating quality threshold for finding best trim lengths. Needed only for Miseq. Defaults to 25.
#' @param max.percent.dropped A positive integer indicating maximal percentage of dropped reads for finding best trim length algorithm. Needed only for MiSeq. Defaults to 20.
#' @param nbases A positive integer indicating how many bases should be used in error learning process. Defaults to 10e8.
#' @param min.length A positive integer. Minimal length of single-end reads. Defaults to 250.
#' @param amplicon.length A positive integer. Length of an analyzed amplicon. Defaults to 450 (V3-V4 16S rRNA).
#' @param min.overlap A positive integer. Minimal length of overlap allowing reads to be merged. Defaults to 50.
#' @param use.mean.quality A logical indicating if mean quality should be used in finding best trim lengths algorithm. Defaults to FALSE.
#' @param reads A string indicating which reads should be processed. One of 'both', 'forward' (R1), 'reverse' (R2). Defaults to 'both'.
#' @param prepareSRA A logical indicating if table for SRA submission should be prepared. Defaults to FALSE.
#' @param SRA.template A path to a TSV template file. If NULL the one in package is used. Defaults to NULL.
#' @param BioSample.template A path to a TSV file with BioSample accessions. Defaults to NULL.
#' @param BioProject.accession A string indicating BioProject accession of a project from which processed samples come. Needed only if prepareSRA=TRUE. Defaults to ''.
#' @param instrument.name A string indicating which sequencing machine was used to generate analyzed reads. Defaults to 'NovaSeq'.
#' @param kit.version A string indicating which version of Illumina's sequencing kit was used to produce analyzed reads. Defaults to 'v1.5'.
#' @param multithread A logical indicating if multiple processors should be used or a positive integer (a number of CPUs to use). Defaults to TRUE (all processors used).
#' @param huge A logical indicating that samples should be processed separately. Usually needed when processing a full NovaSeq run. Defaults to TRUE.
#'
#' @returns List with two elements: seqtab and report.
#'
#' @importFrom dada2 learnErrors filterAndTrim derepFastq dada mergePairs makeSequenceTable removeBimeraDenovo
#' @export denoiseFastQSet

denoiseFastQSet <- function(directory=".",
                            novaseqLoess=NULL,
                            save = TRUE,
                            basename = 'bacteria',
                            pattern = "_R1_001.fastq.gz", # pattern for R1 files, one for R2 will be generated
                            filtered = F, # if TRUE results from previously performed filtering are used
                            error.rates = NULL, # rds object with pre-learned error rates
                            sample.name.sep = "_", # a character on which to split R1 reads filenames to obtain sample names, samplenames will be substrings up to this pattern
                            truncLen = 'auto', # truncation length for trimming, can be 'auto', for NovaSeq generally no trimming is required (either c(0,0) or c(250,250))
                            pool = F, # for more rare taxa set to 'T' or 'pseudo'
                            maxEE = c(2,2),
                            maxN = 0,
                            quality.threshold = 25,
                            max.percent.dropped = 20,
                            nbases = 1e8, # number of bases for learning error rates
                            min.length = 250, # if single-end reads are used
                            amplicon.length = 450,
                            min.overlap = 50,
                            use.mean.quality = F,
                            reads = "both",
                            prepareSRA = F,
                            SRA.template = NULL, # path to a TSV template file, if NULL the one from 'data' will be used
                            BioSample.template = NULL, # name of MIMARKS package
                            BioProject.accession = "",
                            instrument.name = "NovaSeq",
                            kit.version = "1.5",
                            multithread = T,
                            huge = F # if set to 'T' samples will be processed one by one, so that a derep object would not eat too much memory
){


  if(!reads %in% c('both', 'forward', 'reverse')){
    warning("Reads must be one of 'both', 'forward' or 'reverse', unknown value encountered, assuming paired reads ('both')")
    reads = "both"
  }


  # Specify the full path to the fnFs and fnRs
  path <- normalizePath(directory)
  filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory

  if(reads %in% c("both", "forward")){
    fnFs <- sort(list.files(directory, pattern=pattern)) # list of R1 files
    sample.names <- sapply(strsplit(fnFs, sample.name.sep), `[`, 1) # splitting on sample.name.sep and taking the first part
    fnFs <- file.path(path, fnFs)
    filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
    names(filtFs) <- sample.names
  }else{
    fnFs <- NULL
    filtFs <- NULL
  }

  if(reads %in% c("both", "reverse")){
      patternR <- sub("1", "2", pattern) # first '1' changed to '2'
      fnRs <- sort(list.files(directory, pattern=patternR)) # list of R2 files
      sample.names <- sapply(strsplit(fnRs, sample.name.sep), `[`, 1) # splitting on sample.name.sep and taking the first part
      fnRs <- file.path(path, fnRs)
      filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
      names(filtRs) <- sample.names
  }else{
    patternR <- NULL
    fnRs <- NULL
    filtRs <- NULL
  }

  if(instrument.name == "NovaSeq"){
    truncLen <- c(0, 0)
  }else if(instrument.name == 'MiSeq'){
    if(length(truncLen) == 1 && truncLen == 'auto'){
      truncLen <- mgmisc1::find_best_trim_lengths(directory=directory,
                                         forward.pattern=pattern,
                                         reverse.pattern=patternR,
                                         quality.threshold = quality.threshold,
                                         max.percent.dropped = max.percent.dropped,
                                         use.mean = use.mean.quality,
                                         amplicon.length = amplicon.length,
                                         min.overlap = min.overlap,
                                         min.length = min.length)
    }
  }

  fname <- paste0(basename, "_trimming_report.csv")
  if(!filtered){
    message("denoiseFastqSet: filtering raw reads")
    out <- dada2::filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=maxN, maxEE=maxEE, truncLen=truncLen, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=multithread)
    out <- as.data.frame(out)
    colnames(out) <- c('input', 'trimmed')
    out$sample.names <- sapply(strsplit(rownames(out), "_"), '[', 1)
    write.table(out, fname, sep="\t")
    message("denoiseFastQSet: filtering done")
  }else{
    message("denoiseFastqSet: using previously filtered reads")
    out <- read.table(fname, header=T, sep="\t")
    out$sample.names <- sapply(strsplit(rownames(out), "_"), '[', 1)
  }


  # remove samples with zero reads
  filtFs <- filtFs[ names(filtFs) %in% out[out$trimmed > 0, ]$sample.names ]
  filtRs <- filtRs[ names(filtRs) %in% out[out$trimmed > 0, ]$sample.names ]
  sample.names <- out[ out$trimmed > 0, ]$sample.names
  message("denoiseFastQSet: ", length(sample.names), " samples left after removing those with no reads")

  fname <- paste0(basename, "_errorrates.rds")
  if(is.null(error.rates)){
    if(!is.null(novaseqLoess)){
      if(!is.null(fnFs)){
        message("denoiseFastqSet: learning error rates for forward reads")
        errF <- dada2::learnErrors(filtFs, multithread=multithread, nbases=nbases, errorEstimationFunction=novaseqLoess, verbose=T)
      }else{
        errF <- NULL
      }
      if(!is.null(fnRs)){
        message("denoiseFastqSet: learning error rates for reverse reads")
        errR <- dada2::learnErrors(filtRs, multithread=multithread, nbases=nbases, errorEstimationFunction=novaseqLoess, verbose=T)
      }else{
        errR <- NULL
      }
    }else{
      if(!is.null(fnFs)){
        message("denoiseFastqSet: learning error rates for forward reads")
        errF <- dada2::learnErrors(filtFs, multithread=multithread, nbases=nbases, verbose=T)
      }else{
        errF <- NULL
      }
      if(!is.null(fnRs)){
        message("denoiseFastqSet: learning error rates for reverse reads")
        errR <- dada2::learnErrors(filtRs, multithread=multithread, nbases=nbases, verbose=T)
      }else{
        errR <- NULL
      }
    }
    e <- list(errF=errF, errR=errR)
    saveRDS(e, fname)
  }else{
    message("denoiseFastqSet: using pre-computed error rates")
    e <- readRDS(fname)
    errF <- e$errF
    errR <- e$errR
  }

  if(huge){
    message("denoiseFastqSet: huge data, processing each sample separately")
    w <- process_separately(filtFs, filtRs, errF, errR, sample.names, reads, multithread)
    seqtab <- w$seqtab
    report <- w$report
  }else{
    message("denoiseFastqSet: small data, processing all samples together")
    w <- process_together(filtFs, filtRs, errF, errR, sample.names, reads, pool, multithread)
    seqtab <- w$seqtab
    report <- w$report
  }

  message("denoiseFastqSet: removing chimeras")
  seqtab.nochim <- dada2::removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=TRUE)

  if(save){
    fname = paste0(basename, ".seqtab.nochim.rds")
    message("denoiseFastq: saving an RDS object with seqtab.nochim in ", fname)
    saveRDS(seqtab.nochim, fname)
  }
  message("denoiseFastQSet: there were ", dim(seqtab.nochim)[1], " samples and ", dim(seqtab.nochim)[2], " non-chimeric ASVs")
  message("denoiseFastQSet: there were ", 100 * sum(seqtab.nochim)/sum(seqtab), "% of non-chimeras")
  report <- cbind(report, rowSums(seqtab), rowSums(seqtab.nochim))
  track <- merge(out, report, by=0, all.x=T )

  if(reads %in% c("forward", "reverse")){
    colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim")
  }else{
    colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
  }

  fname = paste0(basename, ".track.csv")
  write.table(track, fname, sep="\t", col.names=NA)
  outlist <- list(seqtab=seqtab.nochim, track=track)
  return(outlist)
}

process_together <- function(
    filtFs,
    filtRs,
    errF,
    errR,
    sample.names,
    reads = "both",
    pool = F,
    multithread = T){
  if(reads == 'both' | reads == 'forward'){
    message("denoiseFastqSet: dereplicating forward reads")
    derepFs <- dada2::derepFastq(filtFs, verbose=TRUE)
    # Name the derep-class objects by the sample names
    names(derepFs) <- sample.names
    message("denoiseFastqSet: denoising forward reads")
    dadaFs <- dada2::dada(derepFs, err=errF, multithread=multithread, pool=pool)
  }
  if(reads == 'both' | reads == 'reverse'){
    message("denoiseFastqSet: dereplicating reverse reads")
    derepRs <- dada2::derepFastq(filtRs, verbose=TRUE)
    names(derepRs) <- sample.names
    message("denoiseFastqSet: denoising reverse reads")
    dadaRs <- dada2::dada(derepRs, err=errR, multithread=multithread, pool=pool)
  }
  if(reads == 'both'){
    message("denoiseFastqSet: merging forward and reverse reads")
    mergers <- dada2::mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=10, verbose=TRUE)
    message("denoiseFastqSet: making sequence tab from merged reads")
    seqtab <- dada2::makeSequenceTable(mergers)
    report <- data.frame(denoised=sapply(dadaFs, getN), merged=sapply(mergers, getN))
    rownames(report) <- names(dadaFs)
  }else if(reads == 'forward'){
    message("denoiseFastqSet: making sequence table from forward reads")
    seqtab <- dada2::makeSequenceTable(dadaFs)
    report <-  data.frame(denoised=sapply(dadaFs, getN))
    rownames(report) <- names(dadaFs)
  }else{
    message("denoiseFastqSet: making sequence table from reverse reads")
    seqtab <- dada2::makeSequenceTable(dadaRs)
    report <-  data.frame(denoised=sapply(dadaRs, getN))
    rownames(report) <- names(dadaRs)
  }
  outlist <- list(seqtab=seqtab, report=report)
  return(outlist)
}

process_separately <- function(
    filtFs,
    filtRs,
    errF,
    errR,
    sample.names,
    reads,
    multithread){
  report <- NULL
  if(reads == "forward"){
    message("denoiseFastqSet: making seqtab from forward reads only")
    dadaFs <- vector("list", length(sample.names))
    names(dadaFs) <- sample.names
    for(s in sample.names){
      message("denoiseFastqSet: processing forward reads of sample ", s)
      dF <- dada2::derepFastq(filtFs[[s]])
      dadaFs[[s]] <- dada2::dada(dF, err=errF, multithread=multithread)
    }
    seqtab <- makeSequenceTable(dadaFs)
    report <-  as.data.frame(sapply(dadaFs, getN))
    rownames(report) <- names(dadaFs)
  }else if(reads == "reverse"){
    message("denoiseFastqSet: making seqtab from reverse reads only")
    dadaRs <- vector("list", length(sample.names))
    names(dadaRs) <- sample.names
    for(s in sample.names){
      message("denoiseFastqSet: processing reverse reads of sample ", s)
      dR <- dada2::derepFastq(filtRs[[s]])
      dadaRs[[s]] <- dada2::dada(dR, err=errR, multithread=multithread)
    }
    seqtab <- makeSequenceTable(dadaRs)
    report <-  data.frame(denoised=sapply(dadaRs, getN))
    rownames(report) <- names(dadaRs)
  }else if(reads == "both"){
    message("denoiseFastqSet: making seqtab from merged forward and reverse reads")
    mergers <- vector("list", length(sample.names))
    names(mergers) <- sample.names
    dadas <- vector("list", length(sample.names))
    names(dadas) <- sample.names
    for(s in sample.names) {
      cat("denoiseFastqSet: processing sample:", s, "\n")
      derepF <- dada2::derepFastq(filtFs[[s]])
      dadaFs <- dada2::dada(derepF, err=errF, multithread=multithread)
      derepR <- dada2::derepFastq(filtRs[[s]])
      dadaRs <- dada2::dada(derepR, err=errR, multithread=multithread)
      merger <- dada2::mergePairs(dadaFs, derepF, dadaRs, derepR)
      dadas[[s]] <- dadaFs
      mergers[[s]] <- merger
    }
    rm(derepF); rm(derepR)
    # Construct sequence table
    seqtab <- dada2::makeSequenceTable(mergers)
    report <-  data.frame(denoised=sapply(dadas, getN), merged=sapply(mergers, getN))
    rownames(report) <- names(dadaFs)
  }
  outlist <- list(seqtab=seqtab, report=report)
  return(outlist)
}


getN <- function(x) sum(dada2::getUniques(x))
