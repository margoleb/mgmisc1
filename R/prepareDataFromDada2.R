#' prepareDataFromDada2
#'
#' Prepare an mg-class object based on dada2 results.
#'
#' @param seqtab A named integer matrix with samples in rows and ASVs in columns. Column names should be actual ASVs' sequences. Chimeras should be removed prior to use of prepareDataFromDada2. Such a matrix can be produced with dada2. Defaults to seqtab.nochim.
#' @param taxa A data frame with ASVs in rows and taxlevels i columns or a path to an RDS file containing such an object. Can be produced with dada2. Defaults to NULL.
#' @param tree A phylo-class object or a path to an RDS or nexus file containing one. Defaults to NULL.
#' @param classification.db A path to database allowing classification with dada2's assignTaxonomy(). Defaults to NULL.
#' @param design A data frame with samples as rows and variables as columns or a path to a file containing such table.
#' @param design.sep A string indicating separator needed to read a design file. Defaults to '\\t'.
#' @param design.rowname A string indicating which column in design should be used as row names. Defaults to NULL.
#' @param remove.rare An integer indicating minimal acceptable abundance of kept ASVs. Defaults to 2 (i.e. singletons and doubletons are dropped). When working with high coverage data (e.g. off of a NovaSeq) it should be increased. At 100 000 reads/sample 5 is ok, and for 1 M seqs/sample 10 should be used.
#' @param remove.taxa A vector of strings with contaminating taxa names. Defaults to c('Chloroplast', 'Mitochondria').
#' @param keep.only A string indicating which taxonomic group should be kept. Defaults to 'Bacteria'.
#' @param generate.tree A logical indicating if a phylogentic tree should be computed. Defaults to TRUE.
#' @param alignment.algorithm A string indicating how an alignment should be computed. Current possibilities: 'mafft', 'mothur', 'muscle'. If a reference alignment is used it should be set to 'mothur'. Defaults to 'mothur'.
#' @param tree.algorithm A string idicating how a tree should be computed. Current possibilities: 'clearcut', 'fasttree', 'raxml'. Defaults to 'fasttree'.
#' @param pairwise.alignment A logical indicating if pairwise alignment should be used for tree generation. Doesn't make sense with reference alignment. Defaults to FALSE.
#' @param reference.alignment A path to a reference alignment, such as SILVA. MOTHUR must be used to compute alignment. Defaults to ''.
#' @param start NULL or an integer indicating beginning of an alignment region to keep. Makes sense only for alignments to a reference. If NULL the value is determined so that at least 90\% of sequences are kept. Defaults to NULL.
#' @param stop NULL or an integer indicating end of an alignment region to keep. Makes sense only for alignments to a reference. If NULL the value is determined so that at least 90\% of sequences are kept. Defaults to NULL.
#' @param rarefy Number of rarefaction iterations to average over. If 0 no rarefying is performed. Defaults to 100.
#' @param depth 'auto' or an integer indicating rarefaction depth (i.e. how many randomly picked sequences from each sample to keep). If 'auto' depth is calculated so that information content is maximized. Defaults to 'auto'.
#' @param sample.weight An integer indicating how many sequences is 'worth' one sample, in other words, how much the algorithm should try to keep samples at the cost of decreasing depth. Defaults to 3000.
#' @param min.depth An integer indicating minimal acceptable rarefaction depth. Defaults to 500.
#' @param step An integer indicating interval at which information content will be assessed when optimizing rarefaction depth. Defaults to 10.
#' @param max.samples.discarded A integer indicating maximal percent of samples eliminated due to too low depth. Defaults to 10 (meaning that at least 90\% of samples must be kept).
#' @param basename A string being the common part of file names. Defaults to ''.
#' @param rand.seed An odd integer being a random numbers generator seed. Defaults to 667 (a new better beast).
#' @param fasta.file A path to a fasta file with ASVs sequences. If NULL the file will be generated. Defaults to NULL.
#' @param count.file A path to a count_table file with abundances of ASVs. If NULL the file will be generated. Defaults to NULL.
#' @param OTU A logical or a vector of numbers greater than 0 and lower than 1. If TRUE OTUs will be generated for 0.03 dissimilarity, if a vector of numbers is given OTUs will be computed for dissimilarities in this vector. Defaults to FALSE.
#' @param generate.picrust A logical indicating if PICRUSt2 analysis should be performed. Defaults to TRUE.
#' @param picrust.filename.base A string which will be a common part of names of PICRUSt2-produced files. Defaults to 'picrust'.
#' @param conda A path to conda. Defaults to 'conda'.
#' @param condaenv A string indicating name of a conda environment in which PICRUSt2 is installed. Defaults to 'picrust'.
#' @param picrust A path to PICRUSt2 script. Defaults to 'picrust2\_pipeline.py'.
#' @param picrust.outdir A path to a directory in which PICRUSt2 results will be stored.
#' @param plot.nsti A logical indicating if a plot of NSTI should be generated.
#' @param plot.nsti.var A string indicating which variable should be used to plot NSTI. Defaults to NULL.
#' @param marker.seq A string indicating what kind of marker sequence the experiment uses. Defaults to '16S'.
#' @param save A logical indicating of the resulting object should be saved as an RDS file.
#' @param mothur A path to MOTHUR executable. Defaults to 'mothur'.
#' @param mafft A path to MAFFT executable. Defaults to 'mafft'.
#' @param muscle A path to MUSCLE executable. Defaults to 'muscle'.
#' @param fasttree A path to FastTree executable. Defaults to 'FastTreeMP'.
#' @param raxml A path to RAxML executable. Defaults to 'raxml'.
#' @param description A string with description of the analyzed experiment.
#' @param processors An integer indicating number of processors to use. Defaults to 1.
#' @param ... Other arguments passed on to worker functions.
#'
#' @returns An mg-class object.
#' @export prepareDataFromDada2
#' @import dada2 phangorn vegan GUniFrac

prepareDataFromDada2 <- function(
    seqtab = seqtab.nochim, # either an object or a path to an RDS file
    taxa = NULL, # either an object or a path to an RDS file
    tree = NULL, # either an object or a path to an RDS or nexus file
    classification.db = NULL,
    design = NULL,
    design.sep = "\t", # separator for reading in a design file
    design.rowname = NULL, # which column of a desing table should be used as rownames, needs to contain unique values
    remove.rare = 2, # removing global singletons and doubletons
    remove.taxa = c("Chloroplast", "Mitochondria"),
    keep.only = "Bacteria",
    generate.tree = TRUE,
    alignment.algorithm = "muscle",
    tree.algorithm = "fasttree", #
    pairwise.alignment = FALSE,
    reference.alignment="",
    start = NULL,
    stop = NULL,
    rarefy = 100, # number of rarefaction iterations, if 0 or FALSE no rarefaction will be made
    depth = 'auto', # rarefaction depth (a positive integer) or 'auto'
    sample.weight=3000, # how many times a sample is more important than a sequence/read
    min.depth=500,
    step=10, # interval at which rarefaction calculations are performed
    max.samples.discarded=10, # maximal percentage of the number of samples discarded
    basename = "bacteria",
    rand.seed = 667,
    fasta.file = "seqs.fasta", # name of a fasta file with ASVs sequences, ALWAYS NOT ALIGNED
    count.file = "seqs.count_table", # name of a count file with ASVs number per sample
    OTU = FALSE,
    generate.picrust = FALSE, # generate PICRUSt2 predictions?
    picrust.filename.base = "picrust",
    conda = "conda", # path to conda
    condaenv = "picrust2", # name of conda environment in which PICRUSt2 is installed
    picrust = "picrust2_pipeline.py",
    picrust.outdir = "picrust2_out",
    plot.nsti = FALSE,
    plot.nsti.var = NULL,
    marker.seq = "16S",
    save = TRUE, # save output as an rda object?
    mothur = "mothur",
    mafft = "mafft",
    muscle = "muscle",
    fasttree = "FastTreeMP",
    raxml = "raxml",
    description = "",
    processors = 1,
    directory = ".",
    ...) # other arguments passed to worker functions

{
  directory <- normalizePath(directory)

  outlist <- mgmisc1::new_mg() # an empty mg class object is created

  outlist$metadata$mgmisc.version <- packageVersion("mgmisc1")
  outlist$metadata$mothur.path <- mothur
  outlist$metadata$muscle.path <- muscle
  outlist$metadata$mafft.path <- mafft
  outlist$metadata$fasttree.path <- fasttree
  outlist$metadata$raxml.path <- raxml
  outlist$metadata$conda.path <- conda
  outlist$metadata$picrust.condaenv <- condaenv
  outlist$metadata$picrust.path <- picrust

  if(is.character(design)){
    design <- read.table(design, header=T, sep=design.sep, dec=".")
  }

  if(generate.tree | !OTU){

  }

  if(is.null(design.rowname)){
    message("prepareDataFromDada2: no name of column in the design was given as samplenames (they must match rownames of the seqtab)")
    stop()
  }else if(! design.rowname %in% colnames(design)){
    message(paste0("prepareDataFromDada2: ", design.rowname, " is not a name of a column in design"))
    stop()
  }else if(length(rownames(seqtab) %in% design[[design.rowname]]) == 0 ){
    message("prepareDataFromDada2: samplenames in the design do not match row names of the seqtab")
  }
  if( !is.numeric(remove.rare) & remove.rare != FALSE){
    message("prepareDataFromDada2: remove.rare must be either FALSE or a non-negative integer")
    stop()
  }
  if( !is.numeric(rarefy) & rarefy != FALSE ){
    message("prepareDataFromDada2: number of rarefaction iterations either must be a positive integer or 'rarefy' must be FALSE")
    stop()
  } else if( !is.numeric(depth) ){
    if( depth != 'auto'){
      message("prepareDataFromDada2: rarefying depth must be either 'auto' or an integer greater than 50")
      stop()
    }else{
      depth <- find_best_rarefaction_depth(featuretab=seqtab, sample.weight=sample.weight, min.depth=min.depth, step=step, max.samples.discarded=max.samples.discarded)
      message("prepareDataFromDada2: calculated best rarefaction depth: ", depth)
    }
  }else{
    if(depth < 50){
      message("prepareDataFromDada2: rarefying depth must be an integer greater than 50")
      stop()
    }
  }
  outlist$metadata$rarefy <- rarefy
  outlist$metadata$depth <- depth

  if( !is.numeric(rand.seed) & rand.seed %% 2 != 1){
    stop("prepareDataFromDada2: random numbers generator seed must be an odd integer")
  }
  if(is.null(fasta.file)){
    fasta.file <- paste0(basename, ".fasta")
    fasta.file <- file.path(directory, fasta.file)
  }
  if(is.null(count.file)){
    count.file <- paste0(basename, ".count_table")
    count.file <- file.path(directory, count.file)
  }

  if( !is.null(tree) & class(tree) == 'phylo'){
    message("prepareDataFromDada2: not generating tree, as it was given")
    generate.tree <- FALSE
    tree.file <- paste0(basename, ".tre")
    tree.file <- file.path(directory, tree.file)
    message("preparingDataFromDada2: saving tree in ", tree.file)
    ape::write.tree(tree, tree.file)
  }

  if(reference.alignment == ""){
    message("prepareDataFromDada2: no reference alignment given, pairwise alignment will be used")
    pairwise.alignment <- TRUE
    outlist$metadata$pairwise.alignment <- TRUE
  }

  set.seed(rand.seed)
  outlist$metadata$creation.date <- date()
  outlist$metadata$marker.seq <- marker.seq
  outlist$metadata$rand.seed <- rand.seed
  outlist$metadata$working.dir <- directory
  outlist$metadata$classification.db <- normalizePath(classification.db)
  outlist$metadata$description <- description
  outlist$metadata$basename <- basename
  outlist$metadata$reference.alignment <- normalizePath(reference.alignment)
  outlist$metadata$tree.algorithm <- tree.algorithm
  outlist$metadata$generate.tree <- generate.tree
  outlist$metadata$remove.rare = remove.rare
  outlist$metadata$remove.taxa = remove.taxa
  outlist$metadata$keep.only = keep.only



  if( !is.null(design) ){
    if(is.character(design)){
      message("prepareDataFromDada2: reading design file ", design)
      design <- read.table(design, header=T, sep=design.sep, dec=".")
    }
    rownames(design) <- design[[design.rowname]]
    design <- design[ rownames(design) %in% rownames(seqtab), ]
    seqtab <- seqtab[ rownames(seqtab) %in% rownames(design), ]
    if(nrow(design) == 0){
      stop("prepareDataFromDada2: there were no matching sample names in seqtab and design")
    }
    design <-design[ order(rownames(design)), ]
    seqtab <- seqtab[ order(rownames(seqtab)), ]
    stopifnot(identical(rownames(seqtab), rownames(design))) # must be TRUE!
    if(is.null(design$samplename)){
      design$samplename <- rownames(design)
    }
  }else{
    design <- as.data.frame(rownames(seqtab))
    colnames(design) <- c("samplename")
    rownames(design) <- rownames(seqtab)
  }

  outlist$original.data$sampledata <- design
  outlist$original.data$seqtab <- seqtab

  if( is.null(taxa) & (!is.null(remove.taxa) | !is.null(keep.only)) ){
    if(!is.null(classification.db)){
      message("prepareDataFromDada2: classifying sequences from original seqtab")
      taxa <- dada2::assignTaxonomy(seqtab, classification.db, multithread=T)
      if(save){
        fname <- paste0(basename, ".taxa.rds")
        saveRDS(taxa, fname)
      }
    } else {
      message("prepareDataFromDada2: no object with taxonomy was given, and the sequences cannot be classfied as classification.db is NULL")
    }
  }

  outlist$original.data$featureannot <- taxa

  message("prepareDataFromDada2: generating fasta and count files")
  f <- mgmisc1::generateFasta(seqtab, fasta.file=fasta.file, count.file=count.file) # generated files are named 'seqs.fasta' and 'seqs.count_table'
  outlist$original.data$fasta.file <- f$fasta.file
  outlist$original.data$fasta <- f$fasta
  outlist$original.data$count.file <- f$count.file
  outlist$original.data$count_table <- f$count_table
  message("prepareDataFromDada2: original fasta in ", outlist$original.data$fasta.file, " original count table in ", outlist$original.data$count.file)

  featuretab <- as.data.frame(seqtab)
  colnames(featuretab) <- paste("ASV", seq(1, dim(seqtab)[2], 1), sep="")

  featureannot <- as.data.frame(taxa)
  rownames(featureannot) <- paste("ASV", seq(1, dim(seqtab)[2], 1), sep="")

  message("prepareDataFromDada2: cleaning ASV data")
  out <- mgmisc1::cleanData(featuretab=featuretab, featureannot=featureannot, remove.rare=remove.rare, remove.taxa=remove.taxa, keep.only=keep.only)
  outlist$ASVs$featuretab <- out$featuretab
  outlist$ASVs$featureannot <- out$featureannot
  out <- mgmisc1::makeConsistent(featuretab=outlist$ASVs$featuretab,
                                 sampledata=outlist$original.data$sampledata,
                                 featureannot=outlist$ASVs$featureannot,
                                 fasta=outlist$original.data$fasta,
                                 original.fasta.file=outlist$original.data$fasta.file,
                                 new.fasta.file=outlist$original.data$fasta.file,
                                 count_table=outlist$original.data$count_table,
                                 count.file=outlist$original.data$count.file,
                                 write=T)
  outlist$ASVs$featuretab <- out$featuretab
  outlist$ASVs$sampledata <- out$sampledata
  outlist$ASVs$featureannot <- out$featureannot
  outlist$ASVs$fasta <- out$fasta
  outlist$ASVs$fasta.file <- out$fasta.file
  outlist$ASVs$count_table <- out$count_table
  outlist$ASVs$count.file <- out$count.file
  message("prepareDataFromDada2: filenames: fasta - ", outlist$ASVs$fasta.file, " count - ", outlist$ASVs$count.file)
  outlist$metadata$reference.alignment <- reference.alignment
  outlist$metadata$pairwise.alignment <- pairwise.alignment
  outlist$metadata$alignment.trimming <- list(start, stop)
  names(outlist$metadata$alignment.trimming) <- c("start", "stop")
  outlist$ASVs$aligned = F

  if( generate.tree ){
    message("prepareDataFromDada2: generating a tree")
    t <- mgmisc1::generateTree(fasta.file=f$fasta.file,
                               count.file=f$count.file,
                               reference.alignment=reference.alignment,
                               alignment.algorithm = alignment.algorithm,
                               algorithm = tree.algorithm,
                               start=start,
                               stop=stop,
                               pairwise.alignment=pairwise.alignment,
                               processors=processors)
    fasta <- outlist$ASVs$fasta # never aligned
    aligned.fasta.file <- sub("\\.fasta", ".clean.aligned.fasta", f$fasta.file) # aligned at this point
    nonaligned.fasta.file <- sub("\\.fasta", ".clean.fasta", f$fasta.file) # non aligned
    message(paste0("prepareDataFromDada2: renaming ", t$fasta.file, " to ", aligned.fasta.file))
    system2("mv", paste0(t$fasta.file, " ", aligned.fasta.file))
    outlist$ASVs$fasta.file <- aligned.fasta.file
    count_table <- outlist$ASVs$count_table
    count.file <- sub("\\.count_table", ".clean.count_table", f$count.file)
    message(paste0("prepareDataFromDada2: renaming ", t$count.file, " to ", count.file))
    system2("mv", paste0(t$count.file, " ", count.file))
    tree.file <- sub("\\.fasta", ".clean.tre", f$fasta.file)
    message(paste0("prepareDataFromDada2: renaming ", t$tree.file, " to ", tree.file))
    system2("mv", paste0(t$tree.file, " ", tree.file))
    out <- mgmisc1::makeConsistent(featuretab=outlist$ASVs$featuretab,
                                   sampledata=outlist$ASVs$sampledata,
                                   featureannot=outlist$ASVs$featureannot,
                                   tree=t$tree, tree.file=tree.file,
                                   fasta=fasta,
                                   original.fasta.file=outlist$ASVs$fasta.file,
                                   new.fasta.file=aligned.fasta.file,
                                   count_table=count_table, count.file=count.file,
                                   write=T)
    outlist$ASVs$featuretab <- out$featuretab # final cleaned, but NOT rarefied ASV table
    outlist$ASVs$sampledata <- out$sampledata
    outlist$ASVs$featureannot <- out$featureannot
    outlist$ASVs$fasta <- out$fasta
    outlist$ASVs$fasta.file <- out$fasta.file
    outlist$ASVs$aligned <- t$aligned
    outlist$ASVs$tree.file <- out$tree.file
    outlist$ASVs$tree <- out$tree
    outlist$ASVs$count_table <- out$count_table
    outlist$ASVs$count.file <- count.file
    outlist$metadata$alignment.trimming$start <- t$start
    outlist$metadata$alignment.trimming$stop <- t$stop
    message(paste0("Filenames in out: ", out$fasta.file, " ", out$count.file, " ", out$tree.file))
    message(paste0("Filenames in outlist: ", outlist$ASVs$fasta.file, " ", outlist$ASVs$count.file, " ", outlist$ASVs$tree.file))
    if(identical(rownames(outlist$ASVs$featureannot), colnames(outlist$ASVs$featuretab))){
      message("prepareDataFromDada2: featureannot and featuretab are consistent")
    }
  }else{
    if(!is.null(tree)){
      message(paste0("Filenames in outlist: ", outlist$ASVs$fasta.file, " ", outlist$ASVs$count.file, " ", outlist$ASVs$tree.file))
      out <- mgmisc1::makeConsistent(featuretab=outlist$ASVs$featuretab,
                                     sampledata=outlist$ASVs$sampledata,
                                     featureannot=outlist$ASVs$featureannot,
                                     tree=tree, tree.file=tree.file,
                                     fasta=outlist$ASVs$fasta,
                                     original.fasta.file=outlist$ASVs$fasta.file,
                                     new.fasta.file=outlist$ASVs$fasta.file,
                                     count_table=outlist$ASVs$count_table, count.file=outlist$ASVs$count.file,
                                     write=T)
      outlist$ASVs$featuretab <- out$featuretab # final cleaned, but NOT rarefied ASV table
      outlist$ASVs$sampledata <- out$sampledata
      outlist$ASVs$featureannot <- out$featureannot
      outlist$ASVs$fasta <- out$fasta
      outlist$ASVs$fasta.file <- out$fasta.file
      outlist$ASVs$aligned <- FALSE
      outlist$ASVs$tree.file <- out$tree.file
      outlist$ASVs$tree <- out$tree
      outlist$ASVs$count_table <- out$count_table
      outlist$ASVs$count.file <- count.file
      outlist$metadata$alignment.trimming$start <- start
      outlist$metadata$alignment.trimming$stop <- stop
      message(paste0("Filenames in out: ", out$fasta.file, " ", out$count.file, " ", out$tree.file))
      message(paste0("Filenames in outlist: ", outlist$ASVs$fasta.file, " ", outlist$ASVs$count.file, " ", outlist$ASVs$tree.file))
      if(identical(rownames(outlist$ASVs$featureannot), colnames(outlist$ASVs$featuretab))){
        message("prepareDataFromDada2: featureannot and featuretab are consistent")
      }
    }
  }



  if(as.logical(rarefy)){
    message(paste0("prepareDataFromDada2: rarefying ASVs to ", depth, " reads ", rarefy, " times"))
    fasta.file <- sub("\\.fasta", ".rarefied.fasta", outlist$ASVs$fasta.file)
    tree.file <- sub("\\.tre", ".rarefied.tre", outlist$ASVs$tree.file)
    count.file <- sub("\\.count_table", ".rarefied.count_table", outlist$ASVs$count.file)
    message(paste0("prepareDataFromDada2: rarefied fasta file: ", fasta.file))
    message(paste0("prepareDataFromDada2: rarefied count file: ", count.file))
    message(paste0("prepareDataFromDada2: rarefied tree file: ", tree.file))
    outlist$ASVs$rarefied.data <- mgmisc1::rarefyData(featuretab=outlist$ASVs$featuretab,
                                                      featureannot=outlist$ASVs$featureannot,
                                                      sampledata=outlist$ASVs$sampledata,
                                                      tree=outlist$ASVs$tree,
                                                      tree.file=tree.file,
                                                      fasta=f$fasta,
                                                      fasta.file=outlist$ASVs$fasta.file,
                                                      count_table=f$count_table,
                                                      count.file=count.file,
                                                      rarefy=rarefy,
                                                      depth=depth)
  }
  makeOTUs = FALSE
  if(length(OTU) == 1){
    if(is.logical(OTU)){
      if(OTU){
        message("prepareDataFromDada2: constructing OTUs for 0.03 dissimilarity threshold")
        OTU = c(0.03)
        makeOTUs = TRUE
      }else{
        message("prepareDataFromDada2: not constructing OTUs")
        makeOTUs = FALSE
      }
    }else if( is.numeric(OTU) ){
      makeOTUs = TRUE
    }else{
      makeOTUs = FALSE
      message("prepareDataFromDada2: no valid dissimilarity levels (numbers < 1), not constructing OTUs")
    }
  }
  if(makeOTUs){
    OTU <- as.character(OTU)
      for( diss in OTU ){
        if( !is.na(diss) & diss < 1 ){
          fasta.file <- outlist$ASVs$fasta.file
          count.file <- outlist$ASVs$count.file
          start <- outlist$metadata$alignment.trimming$start
          stop <- outlist$metadata$alignment.trimming$stop
          aligned <- outlist$ASVs$aligned
          message(paste0("prepareDataFromDada2: generating OTUs at ", diss, " dissimilarity level from ", fasta.file, " and ", count.file))

          if(generate.tree & !is.null(reference.alignment)){
            message(paste0("preparedDataFromDada2: using pre-computed alignment in ", fasta.file))
          }else{
            message(paste0("prepareDataFromDada2: starting from unaligned sequences"))
          }
          otus <- mgmisc1::constructOTUs( fasta.file=fasta.file,
                                          count.file=count.file,
                                          sampledata=outlist$ASVs$sampledata,
                                          taxonomy=outlist$ASVs$featureannot,
                                          basename=basename,
                                          diss=diss,
                                          aligned=aligned,
                                          pairwise.alignment=pairwise.alignment,
                                          reference.alignment=reference.alignment,
                                          alignment.algorithm = alignment.algorithm,
                                          start=start,
                                          stop=stop,
                                          classification.db=classification.db,
                                          mothur=mothur,
                                          muscle = muscle,
                                          mafft = mafft,
                                          tree.algorithm = tree.algorithm,
                                          fasttree = fasttree,
                                          raxml = raxml,
                                          processors=processors)

          diss <- as.character(diss)
          otutable <- otus[[diss]]$featuretab
          featureannot <- otus[[diss]]$featureannot
          tree <- otus[[diss]]$tree
          fasta <- otus[[diss]]$fasta
          fasta.file <- otus[[diss]]$fasta.file
          tree.file <- otus[[diss]]$tree.file
          message("Done")
          outlist$OTUs$dists[[diss]] <-otus[[diss]]
          if(remove.rare){
            message(paste0("Removing ", diss, " OTUs with abundance less than ", remove.rare, " across all samples"))
            otutab.final <- otutable[, colSums(otutable) >= remove.rare]
            message("prepareDataFromDada2: fasta file: ", outlist$OTUs$dists[[diss]]$fasta.file)
            outlist$OTUs$dists[[diss]] <- makeConsistent(featuretab          = otutab.final,
                                                         sampledata          = outlist$ASVs$sampledata,
                                                         featureannot        = outlist$OTUs$dists[[diss]]$featureannot,
                                                         tree                = outlist$OTUs$dists[[diss]]$tree,
                                                         tree.file           = outlist$OTUs$dists[[diss]]$tree.file,
                                                         fasta               = outlist$OTUs$dists[[diss]]$fasta,
                                                         original.fasta.file = outlist$OTUs$dists[[diss]]$fasta.file,
                                                         new.fasta.file      = outlist$OTUs$dists[[diss]]$fasta.file,
                                                         count_table         = outlist$OTUs$dists[[diss]]$count_table,
                                                         count.file          = outlist$OTUs$dists[[diss]]$count.file,
                                                         write               = T)
          }
          assign("cos", otus, envir=parent.frame())

          if(rarefy){
            message(paste0("Rarefying OTUs to ", depth, " reads ", rarefy, " times"))
            tree.file <- sub("\\.tre", ".rarefied.tre", outlist$OTUs$dists[[diss]]$tree.file)
            fasta.file <- outlist$OTUs$dists[[diss]]$fasta.file
            count.file <- sub("\\.count_table", ".rarefied.count_table", outlist$OTUs$dists[[diss]]$count.file)
            message(paste0("prepareDataFromDada2: fasta file: ", fasta.file))
            out <- mgmisc1::rarefyData(featuretab=outlist$OTUs$dists[[diss]]$featuretab,
                                       sampledata=outlist$original.data$sampledata,
                                       featureannot=outlist$OTUs$dists[[diss]]$featureannot,
                                       tree=outlist$OTUs$dists[[diss]]$tree,
                                       tree.file=tree.file,
                                       fasta=outlist$OTUs$dists[[diss]]$fasta,
                                       fasta.file=fasta.file,
                                       count_table=outlist$OTUs$dists[[diss]]$count_table,
                                       count.file=outlist$OTUs$dists[[diss]]$count.file,
                                       depth=depth,
                                       rarefy=rarefy)
            outlist$OTUs$dists[[diss]]$rarefied.data <- out
          }
        }else{
          message("prepareDataFromDada2: ", paste0(diss, " is not a valid dissimilarity cutoff for OTUs construction (must be a number < 1)"))
        }
      }
    }

  if(generate.picrust){
    if(marker.seq != "16S"){
      message("prepareDataFromDada2: can't generate PICRUSt2 predictions for sequences other than 16S rRNA")
    }else{
      if( outlist$ASVs$aligned){
        message("prepareDataFromDada2: generating degapped clean fasta file from ", outlist$ASVs$fasta.file)
        write.table(outlist$ASVs$fasta, "temp.fasta", quote=F, sep="", col.names=F)
        fasta.file <- "temp.fasta"
        message("prepareDataFromDada2: using sequences in temp.fasta")
      }else{
        fasta.file <- outlist$ASVs$fasta.file
      }
      message("prepareDataFromDada2: generating PICRUSt2 functional predictions")
      message("prepareDataFromDada2: PICRUSt2 predictions for non-rarefied ASVs")
      if(processors <= 8){
        picrust.processors <- processors
      }else{
        picrust.processors <- 8
      }
      if(dir.exists(picrust.outdir)){
        message("prepareDataFromDada2: PICRUSt output directory found, trying to read in computation results")
        p <- readPICRUSt2Dir(picrust.outdir)
        if(!p$success){
          message("prepareDataFromDada2: ouput directory for PICRUSt2 does not contain full computation results, deleting")
          unlink(picrust.outdir, recursive=T)
          message("prepareDataFromDada2: computing PICRUSt2 predictions for non-rarefied data")
          p <- mgmisc1::preparePICRUSt2(featuretab=outlist$ASVs$featuretab, sampledata=outlist$ASVs$sampledata, fasta.file=fasta.file, filename.base=picrust.filename.base, conda=conda, condaenv=condaenv, picrust=picrust, picrust.outdir=picrust.outdir, plot.nsti=plot.nsti, variable=plot.nsti.var, processors=picrust.processors)
        }
        message("prepareDataFromData2: PICRUSt2 predictions for non-rarefied data read in")
      }else{
        message("prepareDataFromDada2: computing PICRUSt2 predictions for non-rarefied data")
        p <- mgmisc1::preparePICRUSt2(featuretab=outlist$ASVs$featuretab, sampledata=outlist$ASVs$sampledata, fasta.file=fasta.file, filename.base=picrust.filename.base, conda=conda, condaenv=condaenv, picrust=picrust, picrust.outdir=picrust.outdir, plot.nsti=plot.nsti, variable=plot.nsti.var, processors=picrust.processors)
      }
      outlist$KOs$featuretab <- p$KOs$featuretab
      outlist$KOs$featureannot <- p$KOs$featureannot
      outlist$KOs$sampledata <- outlist$ASVs$sampledata
      outlist$KOs$nsti <- p$nsti
      outlist$pathways$featuretab <- p$pathways$featuretab
      outlist$pathways$featureannot <- p$pathways$featureannot
      outlist$pathways$sampledata <- outlist$ASVs$sampledata
      if(rarefy){
        message("prepareDataFromDada2: PICRUSt2 predictions for rarefied ASVs")
        picrust.filename.base <- paste0(picrust.filename.base, "_rarefied")
        picrust.outdir <- paste0(picrust.outdir, "_rarefied")
        message("prepareDataFromDada2: generating degapped clean fasta file from ", outlist$ASVs$rarefied.data$fasta.file)
        if(file.exists("temp.fasta")){
          unlink("temp.fasta")
        }
        write.table(outlist$ASVs$rarefied.data$fasta, "temp.fasta", quote=F, sep="", col.names=F)
        fasta.file <- "temp.fasta"
        message("prepareDataFromDada2: using sequences in ", fasta.file)
        if(dir.exists(picrust.outdir)){
          message("prepareDataFromDada2: PICRUSt output directory found, trying to read in computation results")
          p <- readPICRUSt2Dir(picrust.outdir)
          if(!p$success){
            message("prepareDataFromDada2: ouput directory for PICRUSt2 does not contain full computation results, deleting")
            unlink(picrust.outdir, recursive=T)
            message("prepareDataFromDada2: computing PICRUSt2 predictions for rarefied data")
            p <- mgmisc1::preparePICRUSt2(featuretab=outlist$ASVs$rarefied.data$featuretab,
                                          sampledata=outlist$ASVs$rarefied.data$sampledata,
                                          fasta.file=fasta.file,
                                          filename.base=picrust.filename.base,
                                          conda=conda,
                                          condaenv=condaenv,
                                          picrust=picrust,
                                          picrust.outdir=picrust.outdir,
                                          plot.nsti=plot.nsti,
                                          variable=plot.nsti.var,
                                          processors=picrust.processors)
          }
          message("prepareDataFromData2: PICRUSt2 predictions for rarefied data read in")
        }else{
          message("prepareDataFromDada2: computing PICRUSt2 predictions for rarefied data")
          p <- mgmisc1::preparePICRUSt2(featuretab=outlist$ASVs$rarefied.data$featuretab, sampledata=outlist$ASVs$rarefied.data$sampledata, fasta.file=outlist$ASVs$rarefied.data$fasta.file, filename.base=picrust.filename.base, conda=conda, condaenv=condaenv, picrust=picrust, picrust.outdir=picrust.outdir, plot.nsti=plot.nsti, variable=plot.nsti.var, processors=picrust.processors)
        }
        outlist$KOs$rarefied.data$featuretab <- p$KOs$featuretab
        outlist$KOs$rarefied.data$featureannot <- p$KOs$featureannot
        outlist$KOs$rarefied.data$sampledata <- outlist$ASVs$rarefied.data$sampledata
        outlist$KOs$nsti <- p$nsti
        outlist$pathways$rarefied.data$featuretab <- p$pathways$featuretab
        outlist$pathways$rarefied.data$featureannot <- p$pathways$featureannot
        outlist$pathways$rarefied.data$sampledata <- outlist$ASVs$rarefied.data$sampledata
      }
    }
  }

  class(outlist) <- 'mg'
  if(save){
    rdsfile <- paste0(basename, "_mg.rds")
    saveRDS(outlist, rdsfile)
  }
  return(outlist)

}
