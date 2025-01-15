#' @export get_newfastq
#'


get_newfastq <- function(fq_R1, fq_R2, primerF, primerR, pF_len, pR_len, index_len, max_mismatch) {
  # get R1 and R2 reads and quality
  fq_R1read <- ShortRead::sread(fq_R1)
  fq_R1qual <- Biostrings::quality(fq_R1)
  fq_R2read <- ShortRead::sread(fq_R2)
  fq_R2qual <- Biostrings::quality(fq_R2)
  fq_R1id <- ShortRead::id(fq_R1)
  fq_R2id <- ShortRead::id(fq_R2)

  pF <- Biostrings::DNAString(primerF);
  pR <- Biostrings::DNAString(primerR);

  fs <- as.logical(Biostrings::vcountPattern(pF, ShortRead::narrow(fq_R1read, index_len+1, index_len+pF_len), max.mismatch=max_mismatch))

  rs <- as.logical(Biostrings::vcountPattern(pR, ShortRead::narrow(fq_R1read, index_len+1, index_len+pR_len), max.mismatch=max_mismatch))



  fqR1fs <- ShortRead::narrow(fq_R1read[fs], index_len+pF_len+1, Biostrings::width(fq_R1[fs]))
  fqQ1fs <- ShortRead::narrow(fq_R1qual[fs], index_len+pF_len+1, Biostrings::width(fq_R1[fs]))
  fqR2fs <- ShortRead::narrow(fq_R2read[fs], index_len+pR_len+1, Biostrings::width(fq_R2[fs]))
  fqQ2fs <- ShortRead::narrow(fq_R2qual[fs], index_len+pR_len+1, Biostrings::width(fq_R2[fs]))
  fqI1fs <- fq_R1id[fs]
  fqI2fs <- fq_R2id[fs]
  ind_ffs <- ShortRead::narrow(fq_R1read[fs], 1, index_len)
  ind_rfs <- ShortRead::narrow(fq_R2read[fs], 1, index_len)

  fqR1rs <- ShortRead::narrow(fq_R2read[rs], index_len+pF_len+1, Biostrings::width(fq_R2[rs]))
  fqQ1rs <- ShortRead::narrow(fq_R2qual[rs], index_len+pF_len+1, Biostrings::width(fq_R2[rs]))
  fqR2rs <- ShortRead::narrow(fq_R1read[rs], index_len+pR_len+1, Biostrings::width(fq_R1[rs]))
  fqQ2rs <- ShortRead::narrow(fq_R1qual[rs], index_len+pR_len+1, Biostrings::width(fq_R1[rs]))
  fqI1rs <- fq_R1id[rs]
  fqI2rs <- fq_R2id[rs]
  ind_frs <- ShortRead::narrow(fq_R2read[rs], 1, index_len)
  ind_rrs <- ShortRead::narrow(fq_R1read[rs], 1, index_len)

  fqR1 <- ShortRead::append(fqR1fs, fqR1rs)
  fqR2 <- ShortRead::append(fqR2fs, fqR2rs)
  fqQ1 <- ShortRead::append(fqQ1fs, fqQ1rs)
  fqQ2 <- ShortRead::append(fqQ2fs, fqQ2rs)
  fqI1 <- ShortRead::append(fqI1fs, fqI1rs)
  fqI2 <- ShortRead::append(fqI2fs, fqI2rs)
  ind_f <- ShortRead::append(ind_ffs, ind_frs)
  ind_r <- ShortRead::append(ind_rfs, ind_rrs)

  outlist <- list( fqR1=fqR1, fqR2=fqR2, fqQ1=fqQ1, fqQ2=fqQ2, fqI1=fqI1, fqI2=fqI2, ind_f=ind_f, ind_r=ind_r )
  return(outlist)

}
