#' @export rarefyData

rarefyData <- function(
    experiment = NULL,
    featuretab = NULL,
    sampledata = NULL,
    featureannot = NULL,
    tree = NULL,
    tree.file = NULL,
    fasta = NULL,
    fasta.file = "",
    count_table = NULL,
    count.file = NULL,
    rand.seed = "667",
    rarefy = 1000, # number of rarefactions to average over
    depth = 'auto', # a positive integer (depth) or 'auto'
    sample.weight = 3000,
    min.depth = 500,
    max.samples.discarded = 10,
    step = 10
){

  if(!is.numeric(depth)){
    if(depth == 'auto'){
      depth <- mgmisc1::find_best_rarefaction_depth(featuretab=featuretab,
                                                    sample.weight=sample.weight,
                                                    min.depth=min.depth,
                                                    max.samples.discarded=max.samples.discarded,
                                                    step=step)
    }
  }
  if(!is.null(experiment)){
    featuretab <- featuretab(experiment=experiment, feature="ASVs", rarefied=F)
    featureannot <- featureannot(experiment=experiment, feature="ASVs", rarefied=F)
    sampledata <- sampledata(experiment=experiment, feature="ASVs", rarefied=F)
    tree <- tree(experiment=experiment, feature="ASVs", rarefied=F)
    tree.file <- sub("\\.tre.*", ".rarefied.tre", tree.file(experiment=experiment, feature="ASVs", rarefied=F))
    fasta <- fasta(experiment=experiment, feature="ASVs", rarefied=F)
    fasta.file <- sub("\\.fasta", ".rarefied.fasta", fasta.file(experiment=experiment, feature="ASVs", rarefied=F))
    count_table <- count.table(experiment=experiment, feature="ASVs", rarefied=F)
    count.file <- sub("\\.count_table", ".rarefied.count_table", count.file(experiment=experiment, feature="ASVs", rarefied=F))
    if(!is.null(featuretab(experiment=experiment, feature="ASVs", rarefied=T))){
      message("rarefyData: experiment contains rarefied data, they will be overwritten!")
    }
  }
  set.seed(rand.seed)
  temp <- GUniFrac::Rarefy(otu.tab=featuretab, depth=depth)$otu.tab.rff
  for( i in 1:rarefy ){
    temp <- temp + GUniFrac::Rarefy(otu.tab=featuretab, depth=depth)$otu.tab.rff
    message("rarefyData: ", i)
  }
  featuretab.rarefied <- round( temp/rarefy )

  message("rarefyData: making consistent")
  message("rarefyData: file names: ", fasta.file, " ", count.file, " ", tree.file)
  out <- mgmisc1::makeConsistent(featuretab=featuretab.rarefied,
                                 sampledata=sampledata,
                                 featureannot=featureannot,
                                 tree=tree,
                                 tree.file=tree.file,
                                 fasta=fasta,
                                 original.fasta.file=fasta.file,
                                 new.fasta.file = sub("\\.fasta", ".rarefied.fasta", fasta.file),
                                 count_table=count_table,
                                 count.file=count.file,
                                 write=T)
  out$rarefy <- rarefy
  out$depth <- depth

  return(out)
}
