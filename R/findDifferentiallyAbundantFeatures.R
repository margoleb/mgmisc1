#' @export findDifferentiallyAbundantFeatures

findDifferentiallyAbundantFeatures <- function(
    formula,
    featureSet=NULL,
    featuretab,
    sampledata,
    basename = "diff_features.",
    featureannot=NULL,
    mergeby=0, # by which column the results should be merged with annotation
    q.threshold=0.05,
    log2FC.threshold=1,
    method=c('deseq'), # deseq - DESeq2, aldex - ALDEx2,
    consensus=F,
    plot.diff.features=F,
    num.top.features=10, # how many featues to plot
    write.single=F,
    write.final=T,
    mc.samples=128,
    processors=4){

  outlist <- list()
  outlist$deseq <- NULL
  outlist$aldex <- NULL
  outlist$consensus <- NULL

  if(!is.null(featureSet)){
    featuretab <- featureSet$featuretab
    featureannot <- featureSet$featureannot
    sampledata <- featureSet$sampledata
  }

  var1 <- all.vars(as.formula(formula))[1];
  stopifnot( length(all.vars(as.formula(formula))) == 1, all.vars(as.formula(formula)) %in% colnames(sampledata), length(levels(as.factor(sampledata[[var1]]))) >= 2, identical(rownames(featuretab), rownames(sampledata)) );

  message(paste0("findDifferentiallyAbundantFeatures: finding differentially abundant features"))

  asvtabnzv <- mixOmics::nearZeroVar(featuretab);
  if( nrow(asvtabnzv$Metrics) > 0){
    asvtablv <- featuretab[ , -asvtabnzv$Position ];
  } else {
    asvtablv <- featuretab
  }
  message("findDifferentiallyAbundantFeatures: large variance features selected");

  if(method == 'deseq'){
    message("findDifferentialyAbundantFeatures: using DESeq2")
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(t(asvtablv)+1), colData=sampledata, design=as.formula(formula) );
    message("findDifferentiallyAbundantFeatures: DDS object created");
    contr <- combn( levels(as.factor(sampledata[[all.vars(formula)[1]]])), 2 );
    filenameaccu <- paste0(basename, "_", var1, ".diff.csv");
    print( filenameaccu );
    #
    if( processors > 1 ){
      dds <- try(DESeq2::DESeq(dds, parallel=T, BPPARAM=BiocParallel::MulticoreParam(processors), fitType='mean'));
      if(inherits(dds, "try-error")){
        message("findDifferentiallyAbundantFeatures: estimating dispersions error, using gene-wise dispersion estimates")
        dds <- DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(t(asvtablv)+1), colData=sampledata, design=as.formula(formula) );
        dds <- DESeq2::estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds <- DESeq::nbinomWaldTest(dds)
      }
      accuddsres.df <- NULL
      for( i in 1:ncol(contr) ){
        con <- contr[,i];
        constr <- paste0(con[1], "_", con[2]);
        message(paste0("findDifferentiallyAbundantFeatures: analyzing ", constr))
        filename <- paste0(basename, ".", var1, ".", constr, ".diff.csv");
        message("findDifferentiallyAbundantFeatures: filename: ", filename)
        if( !is.null(featureannot) ){
          row <- rep( paste0(con[1], "_", con[2]), ncol(featureannot)+7);
        } else {
          row <- rep( paste0(con[1], "_", con[2]), 7);
        }
        con <- c(var1, con);
        message("findDifferentiallyAbundantFeatures: generating results for ", constr)
        ddsres <- try(DESeq2::results(dds, contrast=con, parallel=T, BPPARAM=BiocParallel::MulticoreParam(processors)));
        if(inherits(ddsres, "try-error")){
          message("findDifferentiallyAbundantFeatures: no results for ", constr)
        }else{
          message(paste0("findDifferentiallyAbundantFeatures: results for ", constr, " generated"))
          ddsres.df <- as.data.frame(ddsres)[ !is.na(ddsres$padj) & ddsres$padj < q.threshold & abs(ddsres$log2FoldChange) > log2FC.threshold, ];
          if(nrow(ddsres.df) == 0){
            message("findDifferentiallyAbundantFeatures: no differentially abundant features in ", constr)
          }else{
            ddsres.df <- ddsres.df[ order(ddsres.df$log2FoldChange), ];
            if( !is.null(featureannot) ){
              ddsres.df <- merge(ddsres.df, featureannot, by.x=0, by.y=mergeby, all.x=T);
            }
            if(is.null(accuddsres.df)){
              if(nrow(ddsres.df) > 0){
                accuddsres.df <- ddsres.df;
                cnames <- colnames(ddsres.df)
                names(row) <- cnames
                accuddsres.df <- rbind(row, accuddsres.df);
                colnames(accuddsres.df) <- cnames
              }
            } else {
              names(row) <- colnames(accuddsres.df)
              accuddsres.df <- rbind(accuddsres.df, row);
              colnames(accuddsres.df) <- cnames
              accuddsres.df <- rbind(accuddsres.df, ddsres.df);
              colnames(accuddsres.df) <- cnames
            }
            if( write.single ){
              write.table(ddsres.df, filename, sep="\t", col.names=NA);
            }
          }
        }
      }
      outlist$deseq <- accuddsres.df
    } else {
      dds <- DESeq2::DESeq(dds);
      for( i in 1:ncol(contr) ){
        con <- contr[,i];
        constr <- paste0(con[1], "_", con[2]);
        filename <- paste0(basename, ".", var1, ".", constr, ".diff.csv");
        if( !is.null(featureannot) ){
          row <- rep( paste0(con[1], "_", con[2]), ncol(featureannot)+7);
        } else {
          row <- rep( paste0(con[1], "_", con[2]), 7);
        }
        con <- c(var1, con);
        ddsres <- DESeq2::results(dds, contrast=con);
        ddsres.df <- as.data.frame(ddsres)[ !is.na(ddsres$padj) & ddsres$padj < q.threshold & abs(ddsres$log2FoldChange) > log2FC.threshold, ];
        ddsres.df <- ddsres.df[ order(ddsres.df$log2FoldChange), ];
        if( !is.null(featureannot) ){
          ddsres.df <- merge(ddsres.df, featureannot, by.x=0, by.y=mergeby, all.x=T);
        }

        if( i == 1 ){
          accuddsres.df <- ddsres.df;
          cnames <- colnames(ddsres.df)
          names(row) <- cnames
          accuddsres.df <- rbind(row, accuddsres.df);
          colnames(accuddsres.df) <- cnames
        } else {
          names(row) <- colnames(accuddsres.df)
          accuddsres.df <- rbind(accuddsres.df, row);
          colnames(accuddsres.df) <- cnames
          accuddsres.df <- rbind(accuddsres.df, ddsres.df);
          colnames(accuddsres.df) <- cnames
        }
   #     write.table(ddsres.df, filename, sep="\t", col.names=NA);
      }
      outlist$deseq <- accuddsres.df
    }
    if( write.final ){
      write.table(accuddsres.df, filenameaccu, sep="\t", col.names=NA);
    }
  } else if( method == 'aldex'){
    message("findDifferentiallyAbundantFeatures: using ALDEx2")
    clr <- ALDEx2::aldex.clr(t(asvtablv), sampledata[[var1]], mc.samples=mc.samples, useMC=T)
    ald <- ALDEx2::aldex.kw(clr, useMC=T)
    ald.significant <- ald[ ald$kw.eBH < q.threshold | ald$glm.eBH < q.threshold, ]
  } else{
    message("findDifferentiallyAbundantFeatures: unknown method: ", method)
  }

  return(outlist);
}

runDESeq2 <- function(
  featureset = NULL,
  variable = NULL

){
  return(NULL)
}

runALDEx2 <- function(
    featureset = NULL,
    variable = NULL
){
  return(NULL)
}

