#' demultiplexFastQ
#'
#' Demultiplex fastq paired fastq files sequenced by Novogene - pools tagged with Novogene's MIDs and within pools 'inner' MIDs. One pool at a time i.e. R1 and R2 files.
#'
#' @param reads1 A path to a fastq file containig forward (R1) reads.
#' @param reads2 A path to a fastq file containing reverse (R2) reads
#' @param primerF A string - forward primer sequence (the one used in the first round of PCR)
#' @param primerR A string - reverse primer sequence (the one used in the first round of PCR)
#' @param design A tsv file with three columns: samplename, index1 (the one at i5 end), index2 (the one at i7 end)
#' @param outdir A path to a directory where resulting files will be stored. Defaults to '.'.
#' @param max_mismatch An integer indicating how many mismatches are acceptable in index sequences. Defaults to 1.
#' @param processors Number of processors to use. Defaults to 1.
#'
#' @returns Data frame with three columns: R1 file name, R2 file name, number of reads. The table contains as many rows as the number of found index combinations (if no reads were found for a combination it is not inculded in the table).
#' @export demultiplexFastQ

demultiplexFastQ <- function(
    reads1, reads2, # paths to R1 and R2 files
    primerF, primerR, # strings, may contain degenerated bases
    design, # path to a tab-separated file with samplenames and barcode sequences for left and right reads, three first columns need to be named 'samplename', 'index_forward' and 'index_reverse'
    outdir = ".",
    max_mismatch = 1,
    processors=1){

  design <- read.table(design, header=T, sep="\t", dec=".")
  destinations <- gsub("-", "_", as.character(design$samplename)) # no dashes ('-') allowed in sample names

  ## get the fastq to split (raise error if fastq untrimmed/not existing)
  fastq_R1 <- reads1
  fastq_R2 <- reads2

  # register parallel backend
  param <- BiocParallel::MulticoreParam(workers = processors)
  if (!BiocParallel::bpisup(param)) {
    BiocParallel::bpstart(param)
    on.exit(BiocParallel::bpstop(param))
  }
  message("demultiplexFastQ: de-multiplexing sequences in ", reads1, " and ", reads2)
  ## filter and write
  info <-
    BiocParallel::bplapply(seq_along(destinations), function(i) {
      split1 <- file.path(outdir, paste0(destinations[i], "_R1.fastq.gz"))
      split2 <- file.path(outdir, paste0(destinations[i], "_R2.fastq.gz"))
      ## stop if the outfiles already exist (otherwise the output would be appended)
      if (file.exists(split1) | file.exists(split2)) {
        warning("demultiplexFastQ: output files already exist! Overwriting...")
        if (file.exists(split1))
          file.remove(split1)
        if (file.exists(split1))
          file.remove(split2)
      }

      ## split files and save number of kept reads
      kept <- splitFastQNovogene(
        fastq_R1,
        fastq_R2,
        primerF = primerF,
        primerR = primerR,
        outfile_R1 = split1,
        outfile_R2 = split2,
        ind_f = design[ i, 2 ],
        ind_r = design[ i, 3 ],
        max_mismatch = max_mismatch
      )
      dfout <- data.frame(R1 = split1,
                          R2 = split2,
                          kept_reads = kept)
      return(dfout)
    }, BPPARAM = param)

}

