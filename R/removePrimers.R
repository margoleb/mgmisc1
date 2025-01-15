#' @export processFastqDirectory
# processes fastq files in a directory, each is split to files containing only one
# amplicon type, files can be placed in a common dir or each kind goes to its own dir
# splits a pair of fastq files containing different amplicon types (e.g. ITS and 16S) into
# files containing only one amplicon type based on primer sequences (amplicons need to be
# short enough for at least 10 3'-end bases of a primer to be sequenced. Primers will be trimmed
# in produced fastq files
#

processFastqDirectory <- function(
    dir = ".",
    pattern = "R1.fastq.gz", #pattern to identify files with forward reads
    sample.name.sep = "_", # a character on which to split R1 reads filenames to obtain samplenames, samplenames will be substrings up to this pattern
    reads = "both", # process R1, R1 or both?
    out.dir = ".",
    split.by.primers = T,
    split.by.taxonomy = F,
    amplicon.names = c(),
    primerF = c(),
    primerR = c(),
    match.length = 15,
    max.mismatch = 2,
    min.trimmed.length = 50,
    reference.dbs = c(), # vector of paths to reference fasta files (formatted for dada2), they must be given in the same order as amplicon.names
    taxonomic.levels = c(), # vector of strings with taxonomic level and a name to be assigned for each database (separated by a ':'), must be of the form 'Kingdom:Bacteria'
    minBoot = 50,
    multithread = T,
    separate.dirs = F
    ){
  report <- c()

  dir <- normalizePath(dir)
  l <- make_file_lists(dir = dir, pattern = pattern, sample.name.sep = sample.name.sep, reads = reads, amplicon.names = amplicon.names)
  to_process <- l$to_process
  sample.names <- l$sample.names
  out.dir <- normalizePath(out.dir)
  if(!dir.exists(out.dir)){
    dir.create(out.dir)
  }
  if(separate.dirs){
    for(n in amplicon.names){
      dpath <- file.path(out.dir, n)
      if(!dir.exists(dpath)){
        dir.create(dpath)
      }
    }
    dpath <- file.path(out.dir, "unassigned")
    if(!dir.exists(dpath)){
      dir.create(dpath)
    }
  }

    message("processFastqDirectory: file lists prepared")

    track <- NULL
    if(split.by.primers){
      for(s in sample.names){
        if(reads == 'both'){
          f <- to_process[rownames(to_process) == s, ]$fnFs
          r <- to_process[rownames(to_process) == s, ]$fnRs
        }else if( reads == 'forward'){
          f <- to_process[rownames(to_process) == s, ]
          r <- NULL
          primerF <- NULL
        }else if( reads == 'reverse'){
          f <- NULL
          r <- to_process[rownames(to_process) == s, ]
          primerR <- NULL
        }
        setwd(dir)
        message("processFastqDirectory: splitting sample ", s, " files: ", f, ", ", r)
        r <- splitFastqAmpliconsByPrimers(dir = dir,
                                            r1 = f,
                                            r2 = r,
                                            out.dir = out.dir,
                                            separate.dirs = separate.dirs,
                                            amplicon.names = amplicon.names,
                                            primerF = primerF,
                                            primerR = primerR,
                                            match.length = match.length,
                                            max.mismatch = max.mismatch,
                                            min.length = min.trimmed.length)
        res <- c(s, r)
        l <- paste(res, collapse ="\t")
        message("processFastqDirectory: ", l)
        if(is.null(track)){
          track <- res
        }else{
          track <- rbind(track, res)
        }
      }
    }

    if(split.by.taxonomy){
      if(split.by.primers){
        p <- paste0("nonassigned.*", pattern)
        l <- make_file_lists(dir=dir, pattern=p, sample.name.sep=sample.name.sep, reads=reads)
        to_process <- l$to_process
        sample.names <- l$sample.names
      }
      t <- NULL
      message("splitting by taxonomy")
      for(s in sample.names){
        message("sample: ", s)
        if(reads == 'both'){
          f <- to_process[rownames(to_process) == s, ]$fnFs
          r <- to_process[rownames(to_process) == s, ]$fnRs
        }else if( reads == 'forward'){
          f <- to_process[rownames(to_process) == s, ]
          r <- NULL
        }else if( reads == 'reverse'){
          f <- NULL
          r <- to_process[rownames(to_process) == s, ]
        }
        message("processFastqDirectory: splitting sample ", s, " files: ", f, ", ", r)
        r <- splitFastqAmpliconsByTaxonomy(dir = dir,
                                           out.dir = out.dir,
                                           r1 = f,
                                           r2 = r,
                                           amplicon.names = amplicon.names,
                                           reference.dbs = reference.dbs, # vector of paths to reference fasta files (formatted for dada2), they must be given in the same order as amplicon.names
                                           taxonomic.levels = taxonomic.levels, # vector of strings with taxonomic level and a name to be assigned for each database (separated by a ':'), must be of the form 'Kingdom:Bacteria'
                                           minBoot = 50,
                                           multithread = multithread,
                                           separate.dirs = separate.dirs
                                           )
        if(is.null(t)){
          t <- r
        }else{
          t <- rbind(t, r)
        }

      }

      if(split.by.primers){
        track <- cbind(track, t)
      }else{
        rownames(t) <- sample.names
        track <- t
      }
    }

  setwd(out.dir)
  return(track)
}

splitFastqAmpliconsByPrimers <- function(
    dir = ".", # path to a directory where fastq files reside
    out.dir = ".", # it should be a normalized path (in an upstream function 'processFastqDirectory')
    separate.dirs = F, # should filtered files be placed in separate directories for each amplicon type? If TRUE the directories need to be created beforehand, their names must be the same as amplicon names
    r1 = "", # name of a fastq file with R1
    r2 = "", # R2
    sample.name.sep = "_S",
    amplicon.names = c(), # vector of strings, e.g. c("16S", "ITS")
    primerF = NULL, # vector of strings with forward primers e.g. c("CAAGGATTACGCAG", "CCTGACGGACGGGATTAT"), primers need to be given in the same order as amplicon names
    primerR = NULL, # as above, but for reverse primers
    match.length = 10, # how long 3' fragment of primers should be looked for?
    max.mismatch = 2, # how many mismatches are allowed in a hit
    min.length = 50, #' minimal length after trimming primers
    append = F
    ){
  if(!is.null(r1)){
    fq_R1 <- ShortRead::readFastq(dirPath = dir, pattern = r1)
  }
  if(!is.null(r2)){
    fq_R2 <- ShortRead::readFastq(dirPath = dir, pattern = r2)
  }
  # get R1 and R2 reads and quality

  if(append){
    mode = 'a'
  }else{
    mode = 'w'
  }

  pat <- paste0(sample.name.sep, ".*")
  if(!is.null(r1)){
    read_length <- Biostrings::width(sread(fq_R1[1]))
    lname <- sub(pat, "", r1)
    p <- paste0(lname, "_")
    r1_rname <- sub(p, "", r1)
  }
  if(!is.null(r2)){
    read_length <- Biostrings::width(sread(fq_R2[1]))
    lname <- sub(pat, "", r2)
    r2_rname <- sub(p, "", r2)
  }

  out <- c()

  for(n in 1:length(amplicon.names)){
#    message("splitFastqAmpliconsByPrimers: getting ", amplicon.names[n], " amplicon")
    if(separate.dirs){
      dname <- file.path(out.dir, amplicon.names[n])
    }else{
      dname <- out.dir
    }
    setwd(dname)
    # prepare primers' fragments for matching
    pF <- Biostrings::DNAString(substr(x = primerF[n], start = nchar(primerF[n])-match.length, stop = nchar(primerF[n])));
    pF_reversed <- Biostrings::reverseComplement(pF)
    pR <- Biostrings::DNAString(substr(x = primerR[n], start = nchar(primerR[n])-match.length, stop = nchar(primerR[n])));
    pR_reversed <- Biostrings::reverseComplement(pR)

    # generate output filenames
    if(!is.null(r1)){
      r1_filename <- paste0(lname, "_", amplicon.names[n], "_", r1_rname)
    }
    if(!is.null(r2)){
      r2_filename <- paste0(lname, "_", amplicon.names[n], "_", r2_rname)
    }

    # select reads with primers
    if(!is.null(r1)){
      reads.to.select <- Biostrings::vcountPattern(pattern=pR_reversed, subject=ShortRead::sread(fq_R1), max.mismatch=max.mismatch, fixed=F)
      nreads.selected <- sum(reads.to.select)
      names(reads.to.select) <- seq(1,length(reads.to.select))
    }else{
      reads.to.select <- Biostrings::vcountPattern(pattern=pF_reversed, subject=ShortRead::sread(fq_R2), max.mismatch=max.mismatch, fixed=F)
      nreads.selected <- sum(reads.to.select)
      names(reads.to.select) <- seq(1,length(reads.to.select))
    }

    if(nreads.selected > 0){
      if(!is.null(r1)){
        outr1 <- fq_R1[as.numeric(names(reads.to.select[ reads.to.select == 1]))]
    #    message("splitFastqAmpliconsByPrimers: R1 selected")
        # trim primers and sequences beyond
        d <- as.data.frame(Biostrings::vmatchPattern(pattern=pR_reversed, subject = ShortRead::sread(outr1), max.mismatch=max.mismatch, fixed=F))
        outr1_trimmed <- IRanges::narrow(outr1, 1, d$start)
    #    message("splitFastqAmpliconsByPrimers: R1 trimmed")
        # filter out short sequences
    #    message("splitFastqAmpliconsByPrimers: filtering out seqs shorter than ", min.length)
        lengths <- IRanges::width(outr1_trimmed)
        names(lengths) <- seq(1, length(outr1_trimmed))
        longer.than.threshold <- as.numeric(names(lengths[lengths > min.length]))
        outr1_trimmed_long <- outr1_trimmed[longer.than.threshold]
    #   message("splitFastqAmpliconsByPrimers: sequences shorter than ", min.length, " nt were filtered out")
        # select nonassigned reads (without primers)
        fq_R1 <- fq_R1[as.numeric(names(reads.to.select[ reads.to.select == 0]))]
        nreads.selected <- length(outr1_trimmed_long)
        # write filtered fastq files
        ShortRead::writeFastq(outr1_trimmed_long, r1_filename, mode=mode)
#        message("splitFastqAmpliconsByPrimers: wrote ", nreads.selected, " reads to ", r1_filename)
        if(!is.null(r2)){
          outr2 <- fq_R2[as.numeric(names(reads.to.select[ reads.to.select == 1]))]
    #      message("splitFastqAmpliconsByPrimers: paired R2 selected")
    #      message("splitFastqAmpliconsByPrimers: lengths determined")
          outr2_trimmed <- IRanges::narrow(outr2, 1, d$start)
#          message("splitFastqAmpliconsByPrimers: R2 trimmed")
          # filter out short sequences
          outr2_trimmed_long <- outr2_trimmed[longer.than.threshold]
          fq_R2 <- fq_R2[as.numeric(names(reads.to.select[ reads.to.select == 0]))]
          ShortRead::writeFastq(outr2_trimmed_long, r2_filename, mode=mode)
#          message("splitFastqAmpliconsByPrimers: wrote ", nreads.selected, " reads to ", r2_filename)
        }
      }else{
        outr2 <- fq_R2[as.numeric(names(reads.to.select[ reads.to.select == 1]))]
#        message("splitFastqAmpliconsByPrimers: R2 selected")
        outr2_trimmed <- IRanges::narrow(outr2, 1, as.data.frame(Biostrings::vmatchPattern(pattern=pF_reversed, subject = sread(outr2), max.mismatch=max.mismatch, fixed=F))$start)
#        message("splitFastqAmpliconsByPrimers: R2 trimmed")
        # filter out short sequences
        outr2_trimmed_long <- outr2_trimmed[longer.than.threshold]
        fq_R2 <- fq_R2[as.numeric(names(reads.to.select[ reads.to.select == 0]))]
        ShortRead::writeFastq(outr2_trimmed_long, r2_filename, mode=mode)
#        message("splitFastqAmpliconsByPrimers: wrote ", nreads.selected, " reads to ", r2_filename)
      }
      out <- c(out, nreads.selected)
    }else{
      message("splitFastqAmpliconsByPrimers: no reads selected for ", amplicon.names[n])
      out <- c(out, 0)
    }
  }

  # names for non-assigned reads
  r1_filename <- paste0(lname, "_nonassigned_", r1_rname)
  r2_filename <- paste0(lname, "_nonassigned_", r2_rname)

  # write non-assigned reads (never in an 'append' mode, the number of unassigned reads may only be reduced)
  if(separate.dirs){
    dname <- file.path(out.dir, "unassigned")
  }else{
    dname <- out.dir
  }
  setwd(dname)

  n_unassigned <- length(fq_R1)
  if(n_unassigned > 0){
    if(!is.null(r1)){
      ShortRead::writeFastq(fq_R1, r1_filename, mode = 'w')
#      message("splitFastqAmpliconsByPrimers: wrote ", n_unassigned, " reads to ", r1_filename)
      if(!is.null(r2)){
        ShortRead::writeFastq(fq_R2, r2_filename, mode = 'w')
#        message("splitFastqAmpliconsByPrimers: wrote ", n_unassigned, " reads to ", r2_filename)
      }
    }else{
      ShortRead::writeFastq(fq_R2, r2_filename, mode = 'w')
#      message("splitFastqAmpliconsByPrimers: wrote ", n_unassigned, " reads to ", r2_filename)
    }
  }

  out <- c(out, n_unassigned)
  names(out) <- c(amplicon.names, "unnassigned")

  return(out)
}

splitFastqAmpliconsByTaxonomy <- function( # for amplicons too long to be split by primers
  dir = ".",
  out.dir = ".",
  r1 = "",
  r2 = "",
  amplicon.names = c(),
  reference.dbs = c(), # vector of paths to reference fasta files (formatted for dada2), they must be given in the same order as amplicon.names
  taxonomic.levels = c(), # vector of strings with taxonomic level and a name to be assigned for each database (separated by a ':'), must be of the form 'Kingdom:Bacteria'
  minBoot = 50,
  separate.dirs = F,
  append = F,
  multithread = T
){
  message("splitFastqAmpliconsByTaxonomy: processign files ", r1, ", ", r2)
  fq_R1 <- ShortRead::readFastq(dir=dir, pattern=r1)
  fq_R2 <- ShortRead::readFastq(dir=dir, pattern=r2)

  seqs <- as.character(sread(fq_R1))
  nreads <- length(seqs)
  message("splitFastqAmpliconsByTaxonomy: ", nreads, " reads to classify")
  lname <- sub("_.*", "", r1)
  p <- paste0(lname, "_")
  r1_rname <- sub(p, "", r1)
  r2_rname <- sub(p, "", r2)

  if(append){
    mode = 'a'
  }else{
    mode = 'w'
  }

  out <- c()

  for(n in 1:length(amplicon.names)){
    if(separate.dirs){
      dname <- file.path(out.dir, amplicon.names[n])
    }else{
      dname <- out.dir
    }
    setwd(dname)
    message("splitFastqAmpliconsByTaxonomy: getting ", amplicon.names[n], " amplicons, placing results in ", dname)

    aname <- amplicon.names[n]
    db <- reference.dbs[n]
    level <- c(strsplit(taxonomic.levels[n], ':')[[1]][1])
    tname <- strsplit(taxonomic.levels[n], ':')[[1]][2]

    r1_filename <- paste0(lname, "_", amplicon.names[n], "_", r1_rname)
    r2_filename <- paste0(lname, "_", amplicon.names[n], "_", r2_rname)
    nreads <- length(seqs)
    if(nreads > 0){
      cat("splitFastqAmpliconsByTaxonomy: classifying ", nreads, " reads ... ")
      res <- as.data.frame(dada2::assignTaxonomy(seqs=seqs, refFasta=db,  minBoot = minBoot, verbose = T, multithread=multithread))
      message("done")
      message("classified ", nrow(res), " reads")
      rownames(res) <- seq(1:nrow(res))
      reads.to.select <- as.numeric(rownames(res[!is.na(res[[level]]) & res[[level]] == tname, ]))

      if(!is.null(reads.to.select) && length(reads.to.select) > 0){
        nreads.selected <- length(reads.to.select)
        outr1 <- fq_R1[reads.to.select]
        outr2 <- fq_R2[reads.to.select]

        # select nonassigned reads (without primers)
        reads.to.select <- as.numeric(rownames(res[!(!is.na(res[[level]]) & res[[level]] == tname), ]))
          fq_R1 <- fq_R1[reads.to.select]
          fq_R2 <- fq_R2[reads.to.select]
          seqs <- as.character(sread(fq_R1))

          # write filtered fastq files
          ShortRead::writeFastq(outr1, r1_filename, mode=mode)
          ShortRead::writeFastq(outr2, r2_filename, mode=mode)
          message("splitFastqAmpliconsByTaxonomy: wrote ", nreads.selected, " reads to ", r1_filename, " and ", r2_filename)
          out <- c(out, nreads.selected)
      }else{
        message("splitFastqAmpliconsByTaxonomy: no reads selected for ", amplicon.names[n])
        out <- c(out, 0)
      }
    }else{
      message("splitFastqAmpliconsByTaxonomy: no reads selected for ", amplicon.names[n])
      out <- c(out, 0)
    }
  }

  if(length(fq_R1) > 0){
    # names for non-assigned reads
    r1_filename <- paste0(lname, "_nonassigned_", r1_rname)
    r2_filename <- paste0(lname, "_nonassigned_", r2_rname)

    # write non-assigned reads (never in an 'append' mode, the number of unassigned reads may only be reduced)
    message("splitFastqAmpliconsByTaxonomy: writing unassigned fastq files: ", r1_filename, ", ", r2_filename)
    ShortRead::writeFastq(fq_R1, r1_filename, mode = 'w')
    ShortRead::writeFastq(fq_R2, r2_filename, mode = 'w')
  }

  out <- c(out, length(fq_R1))
  names(out) <- c(amplicon.names, "nonassigned")


  return(out)
}

make_file_lists <- function(
    dir = ".",
    pattern = "R1.fastq.gz",
    sample.name.sep = "_", # a character on which to split R1 reads filenames to obtain samplenames, samplenames will be substrings up to this pattern
    reads = "both", # process R1, R1 or both?
    amplicon.names=NULL
){
  fnFs <- NULL
  fnRs <- NULL

  path <- normalizePath(dir)
  if(reads %in% c('forward', 'both')){
    fnFs <- sort(list.files(path=dir, pattern=pattern)) # list of R1 files
    for(a in amplicon.names){
      v <- grep(a, fnFs)
      if(!identical(v, integer(0))){
        fnFs <- fnFs[ -v ]
      }
    }
    sample.names <- unique(sapply(strsplit(fnFs, sample.name.sep), `[`, 1)) # splitting on sample.name.sep and taking the first part
    names(fnFs) <- sample.names
  }

  if(reads %in% c('reverse', 'both')){
    patternR <- gsub("R1", "R2", pattern)
    fnRs <- sort(list.files(path=dir, pattern=patternR)) # list of R2 files
    for(a in amplicon.names){
      v <- grep(a, fnRs)
      if(!identical(v, integer(0))){
        fnRs <- fnRs[-v ]
      }
    }
    sample.names <- unique(sapply(strsplit(fnRs, sample.name.sep), `[`, 1))
    names(fnRs) <- sample.names
  }

  to_process <- as.data.frame(cbind(fnFs, fnRs))
  rownames(to_process) <- sample.names

  l <- list(to_process = to_process, sample.names = sample.names)
  return(l)
}
