#' constructOTUs
#'
#' Cluster ASVs into OTUs. Opticlust algorithm (Schloss and Westcott, 2016?) is used. Tree is computed and classification of represetative sequences is performed if needed.
#'
#' @param seqtab A named integer matrix with samples in rows and ASVs in columns. Column names should be actual ASVs' sequences.
#' @param fasta.file A path to a fasta file containing ASVs' sequences. Needs to be given if there is no seqtab. Defaults to 'seqs.fasta'.
#' @param count.file A path to a count_table file with abundances of ASVs. Needs to be given if there is no seqtab. Defaults to 'seqs.count_table'.
#' @param dist.file A path to a distance matrix file in column or phylip format. Optional. Defaults to NULL.
#' @param diss A vector of numbers between 0 and 1 indicating at which dissimilarities OTUs should be constructed. Defaults to c(0.03).
#' @param aligned A logical indicating if sequences in fasta.file are aligned. Defaults to FALSE.
#' @param reference.alignment A path to a reference alignment in fasta format (e.g. SILVA or Greengenes). Defaults to NULL.
#' @param start An integer indicating beginning of a kept region of an alignment or NULL in which case it will be calculated so that min. 75\% of sequences are kept. Defaults to NULL.
#' @param stop An integer indicating end of a kept region of an alignment or NULL in which case it will be calculated so that min. 75\% of sequences are kept. Defaults to NULL.
#' @param pairwise.alignment A logical indicating if pairwise alignment should be performed. Defaults to FALSE.
#' @param alignment.algorithm A string indicating which algorithm should be used for alignment in case no reference alignment is given, and no pairwise alignment is to be performed. One of 'muscle' or 'mafft'. Defaults to 'muscle'.
#' @param generate.tree A logical indicating if a tree should be constructed for representative sequences. Defaults to TRUE.
#' @param tree.algorithm A string indicating an algorithm for tree construction. One of 'clearcut' (relaxed neighbor joining, Shenneman et al. 2006), 'fasttree' (ML by FastTreeMP) or 'raxml' (ML by RAxML). Defaults to 'fasttree'.
#' @param classification.db A path to a classification database suitable for use with dada2's assignTaxoonmy. Defaults to NULL.
#' @param mothur A path to a MOTHUR executable. Defaults to 'mothur'.
#' @param muscle A path to MUSCLE executable. Defaults to 'muscle'.
#' @param mafft A path to MAFFT executable. Defaults to 'mafft'.
#' @param fasttree A path to FastTreeMP executable. Defaults to 'FastTreeMP'.
#' @param raxml A path to RAxML executable. Defaults to 'RAxML'.
#' @param processors An integer indicating how many processors should be used. Defaults to 1.
#'
#' @returns A list of featuresets, one for each of dissimilarities.
#'
#' @export constructOTUs


constructOTUs <- function(
    experiment = NULL,
    seqtab = NULL,
    taxonomy = NULL,
    fasta.file = "seqs.fasta",
    count.file = "seqs.count_table",
    sampledata = NULL,
    dist.file = NULL, # path to a pre-calculated distance file (in column or phylip format)
    basename = "",
    diss = c("0.03"),
    aligned = F,
    reference.alignment = NULL,
    start=NULL,
    stop=NULL,
    pairwise.alignment = FALSE,
    alignment.algorithm = "muscle",
    classification.db = NULL,
    generate.tree = TRUE,
    tree.algorithm = "fasttree",
    mothur = "mothur", # path to mothur executable
    fasttree = "FastTreeMP",
    raxml = "RAxML",
    muscle = "muscle",
    mafft = "mafft",
    processors = 1){

  diss <- as.character(diss)

  if(!is.null(experiment)){
    fasta.file <- experiment$ASVs$fasta.file
    count.file <- experiment$ASVs$count.file
    sampledata <- experiment$ASVs$sampledata
    basename = experiment$metadata$basename
    aligned = experiment$ASVs$aligned
    reference.alignment = experiment$metadata$reference.alignment
    alignment.algorithm = experiment$metadata$alignment.algorithm
    generate.tree = experiment$metadata$generate.tree
    tree.algorithm = experiment$metadata$tree.algorithm
    start=experiment$metadata$alignment_trimming$start
    stop=experiment$metadata$alignment_trimming$stop
    pairwise.alignment = experiment$metadata$pairwise.alignment
    classification.db = experiment$metadata$classification.db
    mothur = experiment$metadata$mothur.path # path to mothur executable
    muscle = experiment$metadata$muscle.path
    mafft = experiment$metadata$mafft.path
    fasttree = experiment$metadata$fasttree.path
    raxml = experiment$metadata$raxml.path
  }

  if( system2(mothur, "--version") != 0 ){
    stop("constructOTUs: ", mothur, " is not a functional mothur binary")
  }

  if(!xor(is.null(fasta.file), is.null(seqtab)) ){
    message("constructOTUs: either a fasta file and a matching count file or a dada2-produced seqtab must be given")
    stop()
  }

  if(is.null(reference.alignment) & !pairwise.alignment){
    message("constructOTUs: no reference alignment was given, setting pairwise.alignment to TRUE")
    pairwise.alignment <- TRUE
  }

  if(is.null(fasta.file) & !is.null(seqtab )){
    message("constructOTUs: generating fasta and count_table")
    temp <- mgmisc1::generateFasta(seqtab)
    fasta.file <- temp[1]
    count.file <- temp[3]
  }

  if(is.null(tree.algorithm)){
    tree.algorithm <- 'fasttree'
  }

  out <- list()

  if( is.null(dist.file) ){
    df <- sub("\\.fasta", ".pick.good.filter.dist", fasta.file)
    if(!file.exists(df)){
      message("constructOTUs: generating distance matrix")
      res <- generateDistMat(
        fasta.file = fasta.file,
        count.file = count.file,
        aligned=aligned,
        start = start,
        stop = stop,
        format="column",
        pairwise.alignment = pairwise.alignment,
        reference.alignment = reference.alignment,
        algorithm = alignment.algorithm,
        mothur = mothur,
        muscle = muscle,
        mafft = mafft,
        processors=processors)
      dist.file <- res$distance.matrix
      count.file <- res$count.file
      alignment <- res$fasta.file
    }else{
      message("constructOTUs: found a distance matrix file ", df, "; using it")
      dist.file <- df
      alignment <- sub("\\.dist", ".fasta", df) # replace with code checking if the file is present and if not generating the alignment
    }
  }else{
    message("constructOTUs: constructing OTUs based on ", dist.file)
    if(grep("phylip", dist.file)){
      stop("Phylip formatted distance matrix is incompatible with the opticlust algorithm")
    }
  }

  # clustering using OptiMcc
  list.file <- sub("\\.dist", ".opti_mcc.list", dist.file)
  if(!file.exists(list.file)){
    message("constructOTUs: clustering based on distance matrix in ", dist.file, " using OptiMCC algorithm")
    cmd <- paste0("\"#cluster(column=", dist.file, ", count=", count.file, ", processors=", processors, "); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    list.file <- stdout[ grep("\\.opti_mcc\\.list", stdout) ]
  }

  # generating OTUs
  for(cutoff in diss){
    cutoff <- as.character(cutoff)
    message(paste0("constructOTUs: constructing OTUs for ", cutoff, " dissimilarity using ", dist.file, " and ", count.file))

    # generating shared OTU table
    message("constructOTUs: preparing OTU table")
    cmd <- paste0("\"#make.shared(list=", list.file, ", count=", count.file, ", label=", cutoff, "); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    shared.file <- stdout[ grep("\\.opti_mcc\\.shared", stdout) ]
    # reading in shared file
    otutab <- as.data.frame(read.table(shared.file, header=T, sep="\t"))
    rownames(otutab) <- otutab$Group
    otutab$Group <- NULL
    otutab$numOtus <- NULL
    otutab$label <- NULL
    out[[cutoff]]$featuretab <- otutab

    # finding representative sequences
    message("constructOTUs: finding representative sequences")
    cmd <- paste0("\"#get.oturep(list=", list.file, ", column=", dist.file, ", count=", count.file, ", cutoff=", cutoff, ", fasta=", alignment, ", rename=T); quit()\"")
    stdout <- system2(mothur, cmd, stdout=T)
    rep.fasta.file <- stdout[grep(paste0(cutoff, "\\.rep\\.fasta"), stdout)]
    rep.count.file <- stdout[grep(paste0(cutoff, "\\.rep\\.count_table"), stdout)]
    message(paste0("constructOTUs: representative sequences in ", rep.fasta.file, " count table: ", rep.count.file))
    ## changing fasta ids from ASVs to OTUs in rep.fasta file
    f <- ShortRead::readFasta(".", rep.fasta.file)
    f <- ShortRead::ShortRead(sread=ShortRead::sread(f), id=Biostrings::BStringSet(sub("\t.*", "", as.character(ShortRead::id(f)))))
    ShortRead::writeFasta(f, rep.fasta.file, mode='w')
    ids <- paste0(">", as.character(ShortRead::id(f)))
    seqs <- as.character(ShortRead::sread(f))
    out[[cutoff]]$fasta <- data.frame(seq=seqs, row.names=ids)
    ct <- read.table(count.file, header=T, sep="\t")
    out[[cutoff]]$count_table <- ct
    out[[cutoff]]$count.file <- rep.count.file
    out[[cutoff]]$fasta.file <- rep.fasta.file

    # classifying OTUs (consensus taxonomy via MOTHUR's classify.otu)
    message("constructOTUs: classifying OTUs")
    # prepare .taxonomy file with ASVs classification, if needed
    taxonomy_file <- paste0(basename, ".taxonomy")
    if(!file.exists(taxonomy_file)){
      tx <- data.frame(asv=rownames(taxonomy), tax=apply(taxonomy[,1:ncol(taxonomy)], MARGIN = 1, paste, collapse=";"))
      tx$tax <- paste0(tx$tax, ";")
      tx$tax <- gsub("NA;", "", tx$tax)
      write.table(tx, taxonomy_file, quote=F, sep="\t", row.names=F, col.names=F)
    }
    out[[cutoff]]$featureannot <- classify_otus(
      list.file = list.file,
      count.file = count.file,
      taxonomy = taxonomy_file,
      dissimilarity = cutoff,
      mothur = mothur
    )

    # generating tree
    if(generate.tree){
      message(paste0("constructOTUs: generating a tree for OTUs at ", cutoff, " dissimilarity level"))
      if(tree.algorithm == 'clearcut'){
        message(paste0("constructOTUs: generating a tree with clearcut"))
        t <- mgmisc1::generateTree(
          fasta.file = rep.fasta.file,
          count.file = rep.count.file,
          aligned = TRUE,
          pairwise.alignment = pairwise.alignment,
          reference.alignment = reference.alignment,
          alignment.algorithm = alignment.algorithm,
          start = start,
          stop = stop,
          distance.file = "", # path to a mothur-produced distance matrix
          algorithm = "clearcut", # one of 'fasttree', 'raxml' or 'clearcut'
          mothur = mothur,  # path to a mothur executable, default assumes that mothur is in $PATH
          mafft = mafft,
          muscle = muscle,
          fasttree = fasttree,
          raxml = raxml,
          processors = processors
        )
        out[[cutoff]]$tree <- t$tree
        out[[cutoff]]$tree.file <- t$tree.file
      }else{
        message(paste0("constructOTUs: generating a tree with ", tree.algorithm))
        t <- mgmisc1::generateTree(
          fasta.file = rep.fasta.file,
          count.file = rep.count.file,
          aligned = TRUE,
          pairwise.alignment = pairwise.alignment,
          reference.alignment = reference.alignment,
          alignment.algorithm = alignment.algorithm,
          start = start,
          stop = stop,
          distance.file = "", # path to a mothur-produced distance matrix
          algorithm = tree.algorithm, # one of 'fasttree', 'raxml' or 'clearcut'
          mothur = mothur,  # path to a mothur executable, default assumes that mothur is in $PATH
          mafft = mafft,
          muscle = muscle,
          fasttree = fasttree,
          raxml = raxml,
          processors = processors
        )
        out[[cutoff]]$tree <- t$tree
        out[[cutoff]]$tree.file <- t$tree.file
      }
    }else{
      message("constructOTUs: no tree generated")
      out[[cutoff]]$tree <- NULL
      out[[cutoff]]$tree.file <- ""
    }
    v <- "dupa"
    assign( v, out, pos=parent.frame())
    message("constructOTUs: making produced featureset consistent")
    out[[cutoff]] <- makeConsistent(
      featuretab = out[[cutoff]]$featuretab,
      featureannot = out[[cutoff]]$featureannot,
      sampledata = sampledata,
      tree = out[[cutoff]]$tree,
      tree.file = out[[cutoff]]$tree.file,
      fasta = out[[cutoff]]$fasta,
      original.fasta.file = out[[cutoff]]$fasta.file,
      new.fasta.file = out[[cutoff]]$fasta.file,
      count_table = out[[cutoff]]$count_table,
      count.file = out[[cutoff]]$count.file
    )
    if(identical(colnames(out[[cutoff]]$featuretab), colnames(out[[cutoff]]$featureannot))){
      message("constructOTUs: featureset for cutoff of", cutoff, " is consistent")
    }
  }
  return(out)
}


classify_otus <- function(
  list.file = NULL,
  count.file = NULL,
  taxonomy.file = NULL,
  taxlevel.names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
  dissimilarity = cutoff,
  mothur = "mothur"
){
 # classifying OTUs with MOTHUR's classify.otu
 cmd <- paste0("\"#classify.otu(list=", list.file, ", taxonomy=", taxonomy.file, ", count=", count.file, ", label=", dissimilarity, ", probs=F); quit()\"")
 stdout <- system2(mothur, cmd, stdout=T)
 ctax.file <- paste0(sub("\\.list", "", list.file), ".", as.character(dissimilarity), ".cons.taxonomy" )

 # reading in the MOTHUR-produced file and transforming tax strings to columns named as taxlevels
 ntaxlevels <- length(taxlevel.names)
 ctax <- read.table(ctax.file, header=T, sep="\t")
 res <- as.data.frame(do.call(rbind, lapply(strsplit(ctax$Taxonomy, ";"), function(x) { length(x) <- ntaxlevels; return(x) })))
 rownames(res) <- ctax$OTU
 colnames(res) <- taxlevel.names

 return(res)
}
