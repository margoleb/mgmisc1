#' @export print.mg


print.mg <- function( x, ...){
  cat(paste0("\nMetagenome experiment created on ", x$creation.date, " \n"))
  cat(paste0("Marker sequence: ", x$marker.seq, "\n"))
  cat(paste0("Random numbers generator seed used: ", x$rand.seed, " \n"))
  cat(paste0("Directory containing associated files: ", x$working.directory, " \n"))
  cat(paste0("Original seqtab with ", nrow(x$original.data$seqtab), " samples and ", ncol(x$original.data$seqtab), " ASVs"))
  cat(paste0("Sequences were classified with ", x$classification.db, "\n"))
  cat(paste0(ncol(x$original.data$sampledata), " variables in sample data table"))

  if(!is.null(x$ASVs)){
    cat(paste0("ASVs with abundance lower than ", x$ASVs$remove.rare, " were removed\n"))
    cat(paste0("Sequences belonging to ", paste(x$ASVs$remove.taxa), " were removed\n"))
    cat(paste0("Only sequences belonging to ", x$ASVs$keep.only, " were kept\n"))
    cat(paste0(nrow(x$ASVs$featuretab), " samples and ", ncol(x$ASVs$featuretab), " left after clean up\n"))
  }


}
