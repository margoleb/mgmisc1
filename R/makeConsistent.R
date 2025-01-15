#' @export makeConsistent

makeConsistent <- function(
    featuretab = NULL,
    sampledata = NULL,
    featureannot = NULL,
    tree = NULL,
    tree.file = NULL,
    fasta = NULL,
    original.fasta.file = "",
    new.fasta.file = "",
    count_table = NULL,
    count.file = "",
    write = T,
    verbose = TRUE
){
  out <- list()
  if(is.null(featuretab)){ message('makeConsistent: feature table must be present'); stop()}
  if(original.fasta.file != "" && new.fasta.file == ""){ new.fasta.file <- original.fasta.file }

  nsamples_orig <- nrow(featuretab)
  nasv_orig <- ncol(featuretab)

  if( is.null(sampledata) & is.null(featureannot) & is.null(tree) ){message('makeConsistent: no data sampledata, featureannot or tree were given, nothing to do'); stop()}

  if( !is.null(sampledata)){
    featuretab <- featuretab[ rownames(featuretab) %in% rownames(sampledata), ] # only those samples, for which there are sampledata
    featuretab <- featuretab[ , colSums(featuretab) > 0 ] # only features with non-zero counts
    sampledata <- sampledata[ rownames(sampledata) %in% rownames(featuretab), ] # only those samples, for which there are ASVs
    featuretab <- featuretab[ order(rownames(featuretab)), ]
    sampledata <- sampledata[ order(rownames(sampledata)), ]
    if( !identical(rownames(featuretab), rownames(sampledata))){
      stop('makeConsistent: something is wrong, there was no possibility of making ASV table and sampledata rownames identical')
    }
    if( nrow(featuretab) == 0 ){
      stop('makeConsistent: no common row names in ASV table and sampledata, check objects (and their row names)')
    }
    nsamples <- nrow(featuretab)
    nasv <- ncol(featuretab)
   }
    out$sampledata <- sampledata
    nsamples <- nrow(featuretab)

    if( !is.null(featureannot)){
      featuretab <- featuretab[ , colnames(featuretab) %in% rownames(featureannot)] # only those ASVs for which there is featureannot
      featureannot <- featureannot[ rownames(featureannot) %in% colnames(featuretab), ]
      featuretab <- featuretab[ , order(colnames(featuretab))]
      featureannot <- featureannot[ order(rownames(featureannot)), ]
      if(!identical(colnames(featuretab), rownames(featureannot))){
        stop('makeConsistent: something went wrong, sets ASV names and taxa names differ')
      }
      if(ncol(featuretab) == 0){
        stop('makeConsistent: no common names in ASVs table and tax table, check objects (OTUs vs ASVs?)')
      }
      nasv <- ncol(featuretab)
      featuretab <- featuretab[ rowSums(featuretab) > 0, ] # only those samples which have any features
    }


    if(!is.null(tree)){
        tree <- phangorn::midpoint(tree) # midpoint rooting
        tree$edge.length[which(is.na(tree$edge.length))] <- 0 # some edge lenghts can be NA, they are coverted to 0
        featuretab <- featuretab[ , colnames(featuretab) %in% tree$tip.label ] # only ASVs that are in the tree
        featureannot <- featureannot[ rownames(featureannot) %in% colnames(featuretab), ]
        tree <- ape::drop.tip(tree, tree$tip.label[!(tree$tip.label %in% colnames(featuretab))]) # only ASVs that are in the featuretab (it can be rarefied)
        if(ncol(featuretab) == 0){
          stop('makeConsistent: no common names in ASV table and tree, check objects (ASVs vs OTUs?)')
        }
        if(write & !is.null(tree.file)){
          message("makeConsistent: new tree written to ", tree.file)
          ape::write.tree(tree, tree.file)
        }
        nasv <- ncol(featuretab)
    }
    out$featuretab <- featuretab
    out$featureannot <- featureannot
    out$tree <- tree
    out$tree.file <- tree.file

    if(!is.null(fasta)){
      fasta$name <- sub("\n$", "", sub(">", "", rownames(fasta)))
      clean.fasta <- fasta[ fasta$name %in% colnames(featuretab), ]
      clean.fasta$name <- NULL
      out$fasta <- clean.fasta
      if( original.fasta.file != ""){
        if(write){
          message("makeConsistent: new fasta will be written to ", new.fasta.file)
          df <- mgmisc1::readFasta(original.fasta.file)
          df$name <- sub("\n$", "", sub(">", "", rownames(df)))
          clean.fasta <- df[ df$name %in% colnames(featuretab), ]
          clean.fasta$name <- NULL
          write.table(clean.fasta, new.fasta.file, quote=F, sep="", col.names=F)
          out$fasta.file <- new.fasta.file
          message("makeConsistent: new fasta written to ", new.fasta.file)
        }
      }
    }else{
      out$fasta <- NULL
      out$fasta.file <- ""
    }

    if(!is.null(count_table)){
      clean.count_table <- count_table[ count_table$Representative_Sequence %in% colnames(featuretab), ]
      if( count.file != "" ){
        message("makeConsistent: new count table will be written to ", count.file)
        if(write){
          write.table(clean.count_table, count.file, quote=F, sep="\t", row.names=F)
          message("makeConsistent: new count table written to ", count.file)
        }
        out$count_table <- clean.count_table
        out$count.file <- count.file
      }
    }else{
      out$count <- NULL
      out$count.file <- ""
    }
    if(verbose){
      message(paste0("makeConsistent: out of ", nsamples_orig, " samples ", nsamples, " were retained.\nmakeConsistent: out of ", nasv_orig, " features ", nasv, " were retained."))
    }


    return(out)
}


