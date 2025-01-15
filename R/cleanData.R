#' cleanData
#'
#' Remove artifacts (rare variants) and unwanted taxa.
#'
#' @param featuretab A named integer matrix with samples in rows and features (ASVs, OTUs, genes, KO functions etc.) in columns. Can be extracted from an mg-class object with featuretab(). Defaults to NULL.
#' @param featureannot A data frame with features in rows and annotations (usually taxonomy or functions) in columns. Can be extracted from an mg-class object with featureannot(). Defaults to NULL.
#' @param remove.rare An integer indicating maximal total abundance of features removed as rare (thus probably artefactual). Defaults to NULL.
#' @param remove.taxa A vector of strings indicating which taxa are contaminants and should be removed. Defaults to NULL.
#' @param keep.only A string indicating a taxon to keep. Defaults to NULL.
#'
#' @returns List with two elements: featuretab and featureannot.
#' @examples NULL
#' @export cleanData

cleanData <- function(
  featuretab = NULL,
  featureannot = NULL,
  remove.rare = NULL,
  remove.taxa = NULL,
  keep.only = NULL
){
  out <- list()
  initial.nASVs <- ncol(featuretab)
  out$remove.rare <- remove.rare
  if(!is.null(remove.rare) & is.numeric(remove.rare) & remove.rare > 0){
    message(paste0("cleanData: removing ASVs with abundance less than ", remove.rare, " across all samples"))
    featuretab <- featuretab[ , colSums(featuretab) >= remove.rare ]
    featureannot <- featureannot[ rownames(featureannot) %in% colnames(featuretab), ]
    message(paste0("cleanData: after rare taxa removal ", nrow(featureannot), " ASVs remained"))
  } else {
    message("cleanData: not removing rare sequences")
    seqtab.final <- seqtab
  }

  vrs <- colnames(featureannot)
  featureannot$taxonomy <- ""
  for(v in vrs){
    featureannot$taxonomy <- paste0(featureannot$taxonomy, as.character(featureannot[[v]]), sep=";")
  }
  message("cleanData: temporary taxonomy prepared")
  out$removed.taxa <- remove.taxa
  if(!is.null(remove.taxa)){
    pat <- paste(remove.taxa, collapse="|")
    message(paste0("cleanData: removing unwanted taxa: ", pat))
    restaxa <- featureannot[ grep(pat, featureannot$taxonomy, invert=T), ]
    message(paste0("cleanData: after removing unwanted taxa ", nrow(restaxa), " remained"))
    if(nrow(restaxa) == 0){
      stop("cleanData: all taxa removed, something went wrong!")
    }
    featuretab <- featuretab[ , colnames(featuretab) %in% rownames(restaxa)]
  }


  out$kept.taxa <- keep.only
  if(!is.null(keep.only)){
    message(paste0("cleanData: removing taxa which are not ", paste(keep.only)))
    featureannot <- featureannot[grep(keep.only, featureannot$taxonomy),]
    featuretab <- featuretab[ , colnames(featuretab) %in% rownames(featureannot) ]
  }

  featureannot$taxonomy <- NULL

  message(paste0("cleanData: out of ", initial.nASVs, " ", ncol(featuretab), " ASVs left"))

  out$featuretab <- featuretab
  out$featureannot <- featureannot
  return(out)
}
