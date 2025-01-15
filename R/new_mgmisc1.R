#' @export new_mg

new_mg <- function(){
  value <- list()

  creation.date <- date()
  modification.date <- date()
  basename <- ""
  marker.seq <- "16S"
  description <- ""
  rand.seed <- 667
  working.dir <- getwd()
  classification.db <- ""
  mothur <- "mothur"
  mothur.version <- ""
  pairwise.alignment=F
  reference.alignment=""
  conda.path="conda"
  condaenv="picrust2"
  picrust.path="picrust2_pipeline.py"
  PICRUSt2.version=""
  fasttree="FastTreeMP"
  fasttree.version=""
  palette.family="ggsci"
  palette="nrc_npg"
  metadata <- list(creation.date, modification.date, basename, marker.seq, description, rand.seed, working.dir, classification.db, mothur, mothur.version, pairwise.alignment, reference.alignment, conda.path, condaenv, picrust.path, PICRUSt2.version, palette.family, palette)


  seqtab = NULL
  sampledata = NULL
  featureannot = NULL
  tree = NULL
  tree.file = ""
  fasta = NULL
  fasta.file = ""
  count.file = ""
  description=NULL

  original.data <- list(seqtab, sampledata, featureannot, tree, tree.file, fasta, fasta.file, count.file)
  names(original.data) <- c("seqtab", "sampledata", "taxonomy", "tree", "tree.file", "fasta", "fasta.file", "count.file")
  value$original.data <- original.data

  remove.rare = 0
  remove.taxa = c()
  keep.only = ""
  featuretab = NULL
  depth = 0
  rarefy = F
  analyzes = list(alpha.diversity=NULL, beta.diversity=NULL, taxonomy=NULL, diff.features=NULL)
  rarefied <- list(featuretab, sampledata, featureannot, tree, tree.file, fasta, fasta.file, analyzes, depth, rarefy)
  names(rarefied) <- c("featuretab", "sampledata", "featureannot", "tree", "tree.file", "fasta", "fasta.file", "analyzes", "depth", "rarefy")
  ASVs <- list(featuretab, sampledata, featureannot, tree, tree.file, fasta, fasta.file, analyzes, remove.rare, remove.taxa, keep.only, rarefied)
  names(ASVs) <- c("featuretab", "sampledata", "featureannot", "tree", "tree.file", "fasta", "fasta.file", "analyzes", "remove.rare", "remove.taxa", "keep.only", "rarefied.data")
  value$ASVs <- ASVs

  dists <- list()
  OTUs <- list(dists, remove.rare)
  names(OTUs) <- c("dists", "remove.rare")
  value$OTUs <- OTUs

  rarefied.data <- list(featuretab, sampledata, featureannot, fasta, fasta.file, analyzes, depth, rarefy)
  names(rarefied.data) <- c("featuretab", "sampledata", "featureannot", "fasta", "fasta.file", "analyzes", "depth", "rarefy")

  KOs <- list(featuretab, sampledata, description, tree, fasta, fasta.file, analyzes, rarefied.data)
  names(KOs) <- c("featuretab", "sampledata", "featureannot", "tree", "fasta", "fasta.file", "analyzes", "rarefied.data")
  value$KOs <- KOs

  pathways <- list(featuretab, sampledata, description, tree, fasta, fasta.file, analyzes, rarefied.data)
  names(pathways) <- c("featuretab", "sampledata", "featureannot", "tree", "fasta", "fasta.file", "analyzes", "rarefied.data")
  value$pathways <- pathways

  class(value) <- "mg"
  return(value)
}
