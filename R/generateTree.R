#' @export generateTree

generateTree <- function(
    seqtab = NULL,
    fasta.file = "", # path to a fasta file
    count.file = "", # path to a mothur-produced count_table
    aligned = FALSE,
    pairwise.alignment = FALSE,
    reference.alignment = "",
    alignment.algorithm = "muscle",
    start = NULL,
    stop = NULL,
    distance.file = "", # path to a mothur-produced distance matrix
    algorithm = "fasttree", # one of 'fasttree', 'raxml' or 'clearcut'
    mothur = "mothur",  # path to a mothur executable, default assumes that mothur is in $PATH
    mafft = "mafft",
    muscle = "muscle",
    fasttree = "FastTreeMP",
    raxml = "raxml",
    processors = 1) {

  if(is.null(seqtab) & fasta.file == "" & distance.file == ""){
    message("generateTree: either fasta file or seqtab or a distance matrix file must be given")
    stop()
  }

  if(!is.null(seqtab)){
    bname="seqs."
    if(fasta.file == ""){
      message("generateTree: generating fasta and count_table")
      temp <- mgmisc1::generateFasta(seqtab)
      fasta.file <- temp[1]
      count.file <- temp[3]
    }
  }else{
    bname <- gsub("fasta", "", fasta.file)
  }
  result <- list()
  result$start <- start
  result$stop <- stop
  result$dist.phylip.file <- distance.file
  result$mothur <- mothur
  result$mafft <- mafft
  result$muscle <- muscle
  result$raxml <- raxml
  result$fasttree <- fasttree
  result$aligned <- TRUE # resulting fasta file is aligned

    # tree is computed using relaxed neighbor joining algorithm (Sheneman et al.)
    if(algorithm == "clearcut"){
      message("generateTree: computing tree with clearcut")
      if(distance.file == ""){
        if(!aligned & is.null(reference.alignment) & !pairwise.alignment){
          message("generateTree: reference-free alignment will be computed with ", alignment.algorithm)
        }
        tree <- NULL

        message("generateTree: calculating distance matrix")
        dm.res <- mgmisc1::generateDistMat(fasta.file = fasta.file,
                                        count.file = count.file,
                                        start = start,
                                        stop = stop,
                                        aligned = aligned,
                                        pairwise.alignment = pairwise.alignment,
                                        reference.alignment = reference.alignment,
                                        algorithm = alignment.algorithm,
                                        format = "phylip",
                                        mothur = mothur,
                                        mafft = mafft,
                                        muscle = muscle,
                                        processors=processors)
        d <- dm.res$distance.matrix
        result$fasta.file <- dm.res$fasta.file
        result$aligned <- dm.res$aligned
        result$start <- dm.res$start
        result$stop <- dm.res$stop
        result$alignment.algorithm <- alignment.algorithm
        result$mothur.version <- dm.res$mothur.version
        result$mafft.version <- dm.res$mafft.version
        result$muscle.version <- dm.res$muscle.version
      }else{
        message(paste0("generateTree: computing tree from supplied distance matrix file ", distance.file))
        d <- distance.file
        result$aligned <- FALSE
      }
      message(paste0("generateTree: computing tree using distance matrix in ", d))
      cmd <- paste0("\"#clearcut(phylip=", d, "); quit();\"")
      stdout <- system2(mothur, cmd, stdout=T)
      tree.file <- stdout[ grep("Output\ File\ Names", stdout)+1 ]
    }else{
      if(!aligned){
        if(is.null(reference.alignment)){
        message("generateTree: computing alignment using ", alignment.algorithm)
        }else{
          message("generateTree: aligning sequences to ", reference.alignment, " using MOTHUR")
        }
        a <- generateAlignment(
          fasta.file = fasta.file,
          count.file = count.file,
          reference.alignment = reference.alignment,
          algorithm = alignment.algorithm,
          start = start,
          stop = stop,
          mothur = mothur,
          mafft = mafft,
          muscle = muscle,
          processors = processors
        )
        result$start <- a$start
        result$stop <- a$stop
        result$fasta.file <- a$alignment
        result$count.file <- a$count.file
        result$alignment.algorithm <- alignment.algorithm
        result$mothur.version <- a$mothur.version
        result$mafft.version <- a$mafft.version
        result$muscle.version <- a$muscle.version
      }else{
        result$fasta.file <- fasta.file
      }


      if(algorithm == "fasttree"){
        message("generateTree: computing tree using FastTree")
        fasttree_version <- try(system2(fasttree, "-expert", stdout = T))
        if(!inherits(fasttree_version, "try-error")){
          tree.file <- paste0(bname, "fasttree.tre")
          cmd <- paste0(" -fastest -nosupport -nt ", result$fasta.file, " > ", tree.file)
          fasttree_version <- sub("\\ .*", "", sub("Detailed\\ usage\\ for\\ FastTree\\ ", "", fasttree_version[1]))
          stdout <- system2(fasttree, cmd, stdout = T)
        }else{
          stop("generateTree: ", fasttree, " is not a valid FastTreeMP executable")
        }
      }else if(algorithm == "raxml"){
        message("generateTree: computing tree using RAxML")
        raxml_version <- try(system2(raxml, " -v"))
        if(!inherits(raxml_version, "try-error")){
          cmd <- paste0("-m GTRCAT -F -s ", a$alignment, " -n ", tree.file, " -p ", rand.seed, " -T ", processors)
          stdout <- system2(raxml, cmd, stdout = T)
          raxml_version <- sub("This\\ is\\ RAxML\\ version\\ ", "", sub("\\ .*", "", raxml_version[3]))
          tree.file <- paste0("RAxML_result.", tree.file)
        }else{
          stop("generateTree: ", raxml, " is not a valid RAxML executable")
        }
      }else{
        stop("generateTree: unknown tree construction algorithm ", algorithm, ", no tree computed")
      }
    }
    message(paste0("generateTree: tree in ", tree.file))
    t <- ape::read.tree(tree.file)
    message("generateTree: tree generated")
    result$tree <- t
    result$tree.file <- tree.file
    return(result)
}
