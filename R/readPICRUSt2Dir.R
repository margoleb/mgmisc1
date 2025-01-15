#' @export readPICRUSt2Dir


readPICRUSt2Dir <- function(dir=NULL){

  result <- list()
  result$success <- TRUE
  f1 <- paste0(dir, "/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz")
  f2 <- paste0(dir, "/pathways_out/path_abun_unstrat_descrip.tsv.gz")
  f3 <- paste0(dir, "/KO_metagenome_out/weighted_nsti.tsv.gz")

  if(file.exists(f1)){
    picrust_ko <- read.table(f1, sep='\t', dec='.', header=T, row.names=1, quote='"', check.names=F)
    ko_descriptions <- data.frame(row.names=rownames(picrust_ko), fun=rownames(picrust_ko), desc=picrust_ko$description)
    picrust_ko$description <- NULL
    picrust_ko <- t(picrust_ko)
  }else{
    result$success <- FALSE
    return(result)
  }
  # pathways
  if(file.exists(f2)){
    picrust_pathways <- read.table(f2, sep='\t', dec='.', header=T, row.names=1, quote='"', check.names=F)
    pathways_descriptions <- data.frame(row.names=rownames(picrust_pathways), fun=rownames(picrust_pathways), desc=picrust_pathways$description)
    picrust_pathways$description <- NULL
    picrust_pathways <- t(picrust_pathways)
  }else{
    result$success <- FALSE
    return(result)
  }

  # NSTI
  if(file.exists(f3)){
    nsti <- read.table(f3, header=T, sep="\t", dec=".")
  }else{
    result$success <- FALSE
    return(result)
  }


  result$KOs <- list(NULL)
  result$pathways <- list(NULL)
  result$KOs$featuretab <- round(picrust_ko)
  result$KOs$featureannot <- ko_descriptions
  result$pathways$featuretab <- round(picrust_pathways)
  result$pathways$featureannot <- pathways_descriptions
  result$nsti <- nsti

  return(result)
}
