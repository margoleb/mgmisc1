#' @export generateDistMat

# distance matrix calculation for aligned or unaligned fasta sequences, returns filename of a phylip-formatted DM

generateDistMat <- function(
    fasta.file = NULL,
    count.file = NULL,
    aligned = F,
    pairwise.alignment = F,
    reference.alignment = "",
    algorithm = "muscle",
    start = NULL,
    stop = NULL,
    format = "phylip", # either 'phylip' or 'column'
    mothur = "mothur",
    mafft = "mafft",
    muscle = "muscle",
    processors = 1
){
  result <- list()
  result$aligned <- TRUE
  mothur.version <- NULL
  mafft.version <- NULL
  muscle.version <- NULL

  if( aligned ){
    message(paste0("    generateDistMat: calculating distance matrix using aligned sequences"))
    # distance matrix is calculated and saved in the phylip format
    cmd <- paste0("\"#dist.seqs(fasta=", fasta.file, ", output=", format, ", processors=", processors, "); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    d <- stdout[ grep("Output\ File\ Names", stdout)+1 ]
    alignment <- fasta.file
  }else{
    if( pairwise.alignment){
      message("    generateDistMat: calculating distance matrix using pairwise alignment")
      cmd <- paste0(" \"#pairwise.seqs(fasta=", fasta.file, ", output=", format, ", processors=", processors, "); quit();\"")
      stdout <- system2(mothur, cmd, stdout=T)
      d <- stdout[ grep("\\.dist", stdout) ]
      alignment <- NULL
      result$aligned <- FALSE
    }else{
      a <- generateAlignment(
        fasta.file = fasta.file,
        count.file = count.file,
        reference.alignment = reference.alignment,
        algorithm = algorithm,
        start = start,
        stop = stop,
        mothur = mothur,
        mafft = mafft,
        muscle = muscle,
        processors = processors)
      result$start <- a$start
      result$stop <- a$stop
      message("    generateDistMat: computing distance matrix from the final alignment")
      # distance matrix is calculated and saved in the phylip format
      cmd <- paste0("\"#dist.seqs(fasta=", a$alignment, ", output=", format, ", processors=", processors, "); quit()\"")
      stdout <- system2(mothur, cmd, stdout=T)
      d <- stdout[ grep("Output\ File\ Names", stdout)+1 ]
      mothur.version <- a$mothur.version
      mafft.version <- a$mafft.version
      muscle.version <- a$muscle.version
      alignment <- a$alignment
    }
  }
  print(d)
  result$distance.matrix <- d
  result$algorithm <- algorithm
  result$mothur <- mothur
  result$mothur.version <- mothur.version
  result$mafft <- mafft
  result$mafft.version <- mafft.version
  result$muscle <- muscle
  result$muscle.version <- muscle.version
  result$fasta.file <- alignment
  result$count.file <- count.file

  return(result)
}
