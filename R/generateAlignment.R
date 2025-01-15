#' @export generateAlignment

generateAlignment <- function(
    fasta.file = "",
    count.file = "",
    reference.alignment = "",
    start = NULL,
    stop = NULL,
    algorithm = "muscle", # for reference-free alignments: either 'muscle' (muscle -super5) or mafft
    mothur = "mothur",
    mafft = "mafft",
    muscle = "muscle",
    processors = 1){

  outlist <- list()
  outlist$start <- start
  outlist$stop <- stop
  outlist$fasta.file <- fasta.file
  outlist$count.file <- count.file
  outlist$mothur.version <- NULL
  outlist$mafft.version <- NULL
  outlist$muscle.version <- NULL

  if(reference.alignment != ""){
      message("    generateAlignment: calculating distance matrix using reference alignment with Mothur")
    if(system2(mothur, "--version") != 0){
      stop(paste0("prepareDataFromDada2: ", mothur, " is not a valid MOTHUR executable"))
    }else{
      outlist$mothur.path <- mothur
    }
    s <- try(system2(mothur, "--version", stdout=T))
    if(!inherits(s, 'try-error')){
      mothur.version <- strsplit(s[2], "=")[[1]][2]
      outlist$metadata$mothur.version <- mothur.version
    }
      message("    generateAlignment: aligning sequences against ", reference.alignment)
      cmd <- paste0(" \"#align.seqs(fasta=", fasta.file, ", reference=", reference.alignment, ", processors=", processors, "); quit();\"")
      stdout <- system2(mothur, cmd, stdout=T)
      alignment <- stdout[ grep("Output\ File\ Names", stdout)+1 ]
      seqs_to_remove <- stdout[ grep("Output\ File\ Names", stdout)+3 ]
      if(file.exists(seqs_to_remove)){
        message(paste0("    generateAlignment: removing badly aligned seqs listed in ", seqs_to_remove))
        cmd <- paste0(" \"#remove.seqs(fasta=", alignment, ", count=", count.file, ", accnos=", seqs_to_remove, "); quit();\"")
        stdout <- system2(mothur, cmd, stdout=T)
        alignment <- stdout[ grep("Output\ File\ Names", stdout)+1 ]
        count.file <- stdout[ grep("Output\ File\ Names", stdout)+2 ]
        message("    generateDistMat: badly aligned seqs removed")
      }
      # if no desired alignment regions was given (start and stop), it is calculated as 75th percentile of start and stop
      if(is.null(start) | is.null(stop)){
        message("    generateAlignment: determining start and stop positions")
        cmd <- paste0("\"#count.seqs(count=", count.file, "); quit()\"")
        stdout <- system2(mothur, cmd, stdout=T)
        count.file <- stdout[ grep("Output\ File\ Names", stdout)+1 ]
        cmd <- paste0("\"#summary.seqs(fasta=", alignment, ", count=", count.file, ", processors=", processors, "); quit();\"")
        stdout <- system2(mothur, cmd, stdout=T)
        start <- strsplit(stdout[grep("25\\%-tile", stdout)], "\t")[[1]][2]
        stop <- strsplit(stdout[grep("25\\%-tile", stdout)], "\t")[[1]][3]
        message(paste0("    generateAlignment: alignment will be trimmed to ", start, " - ", stop, " position"))
        outlist$start <- start
        outlist$stop <- stop
        outlist$count.file <- count.file
      }
      # screening for sequences covering the desired region of the alignment
      message("    generateAlignment: trimming alignment")
      cmd <- paste0("\"#screen.seqs(fasta=", alignment, ", count=", count.file, ", start=", start, ", end=", stop, ", processors=", processors, "); quit()\"")
      stdout <- system2(mothur, cmd, stdout=T)
      alignment <- stdout[ grep("good.align", stdout) ]
      count.file <- stdout[ grep('good.count_table', stdout)]
      message("    generateAlignment: filtering alignment")
      # all-gap columns and columns containing at least one terminal gap are removed
      cmd <- paste0("\"#filter.seqs(fasta=", alignment, ", vertical=T, trump=., processors=", processors, "); quit()\"")
      stdout <- system2(mothur, cmd, stdout=T)
      alignment <- stdout[ grep("Output\ File\ Names", stdout)+2 ]
      outlist$alignment <- alignment
      outlist$count.file <- count.file
  }else{
    bname <- sub("\\.fasta", "", fasta.file)

    if(algorithm == 'mafft'){
    mafft_version <- try(system2(mafft, "--version", stdout = T))
    if(!inherits(mafft_version, 'try-error')){
      message("generateAlignment: computing a reference-free alignment with MAFFT ", mafft_version)
      alignment <- paste0(bname, "mafft_aln.fasta")
      cmd <- paste0(" --thread ", processors, " fasta.file > ", alignment)
      stdout <- system2(mafft, cmd, stdout = T)
      outlist$mafft.version <- mafft_version
      outlist$mafft.path <- mafft
    }else{
      stop("generateAlignment: ", mafft, " is not a working MAFFT executable")
    }
    }else if(algorithm == 'muscle'){
      muscle_version <- try(system2(muscle, "--version"))
      if(!inherits(muscle_version, "try-error")){
        message("generateAlignment: computing reference-free alignment with MUSCLE ", muscle_version)
        alignment <- paste0(bname, "muscle_aln.fasta")
        outlist$muscle.version <- muscle_version
        outlist$muscle.path <- muscle
        cmd <- paste0(" -threads ", processors, " -super5 ", fasta.file, " -output ", alignment)
        stdout <- system2(muscle, cmd, stdout = T)
      }else{
        stop("generateAlignment: ", muscle, " is not a working MUSCLE executable")
      }
    }else{
      stop("generateAlignment: unknown alignment algorithm ", algorithm)
    }
    outlist$alingment <- alignment
  }

  return(outlist)
}
