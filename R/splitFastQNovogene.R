#' Split paired fastq files according to smaple indexes
#'
#' @param reads1
#' @param reads2
#' @param primerF
#' @param primerR
#' @param outfile_R1
#' @param outfile_R2
#' @param ind_f
#' @param ind_r
#' @param max_mismatch
#'
#' @export splitFastQNovogene
#'
#'
splitFastQNovogene <- function( reads1,
                                reads2,
                                primerF,
                                primerR,
                                outfile_R1,
                                outfile_R2,
                                ind_f,
                                ind_r,
                                max_mismatch ){ # reads1 and reads2 ShortReadQ objects for paired fastq files (order of reads must be the same in both objects) obtained with readFastq of ShortRead

  pF <- toupper(primerF);
  pR <- toupper(primerR);
  indf <- toupper(ind_f);
  indr <- toupper(ind_r);
  indlen <- width(ind_f);

  stream_R1 <- ShortRead::FastqStreamer(reads1)
  stream_R2 <- ShortRead::FastqStreamer(reads2)
  on.exit(close(stream_R1))
  on.exit(close(stream_R2), add = TRUE)

  kept_reads <- 0
  ## filter fastq
  repeat {
    # input chunk
    fq_R1 <- ShortRead::yield(stream_R1, 10000)
    fq_R2 <- ShortRead::yield(stream_R2, 10000)
    if (length(fq_R1) == 0) {
      break
    }
    ## modify R1/R2 as per the protocol
    outlist <- get_newfastq(fq_R1, fq_R2, primerF=pF, primerR=pR, pF_len=width(pF), pR_len=width(pR), index_len=indlen, max_mismatch=max_mismatch)
    ## make a new ShortReadQ object with new header and trimmed reads
    # new fastq R2
    fqid_R2 <- ShortRead::id(fq_R2)
    fqid_R2 <- sub(" ", "_", fqid_R2) # replace space by underscore in readName
    fq_R2new <- ShortRead::ShortReadQ(outlist$fqR2,
                                      outlist$fqQ2,
                                      Biostrings::BStringSet(paste(
                                        outlist$fqI2,
                                        outlist$ind_f,
                                        outlist$ind_r,
                                        sep = "#"
                                      )))

    # new fastq R1
    fqid_R1 <- ShortRead::id(fq_R1)
    fqid_R1 <- sub(" ", "_", fqid_R1) # replace space by underscore in readName

    fq_R1new <- ShortRead::ShortReadQ(outlist$fqR1,
                                      outlist$fqQ1,
                                      Biostrings::BStringSet(paste(
                                        outlist$fqI1,
                                        outlist$ind_f,
                                        outlist$ind_r,
                                        sep = "#"
                                      )))

    # demultiplex using given sample barcodes
    idf_seq <- Biostrings::DNAString(indf)
    idr_seq <- Biostrings::DNAString(indr)
    index_f <- Biostrings::DNAStringSet(outlist$ind_f)
    index_r <- Biostrings::DNAStringSet(outlist$ind_r)
    id2keep <- as.logical( Biostrings::vcountPattern(idf_seq,
                                                     outlist$ind_f,
                                                     max.mismatch = max_mismatch) &
                             Biostrings::vcountPattern(idr_seq,
                                                       outlist$ind_r,
                                                       max.mismatch = max_mismatch) )

    # append to destination
    ShortRead::writeFastq(fq_R1new[id2keep], outfile_R1, "a")
    ShortRead::writeFastq(fq_R2new[id2keep], outfile_R2, "a")
    # add to count
    kept_reads <- kept_reads + sum(id2keep)

  }
  return <- kept_reads;
}
