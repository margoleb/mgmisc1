#' @export plotBetaDiversity
#' @import vegan ggplot2 ggpubr GUniFrac
plotBetaDiversity <- function(
    featureSet=NULL,
    featurename=NULL,
    featuretab=NULL,
    sampledata=NULL,
    featureannot=NULL,
    tree=NULL,
    shape.variable=NULL, # which variable should be mapped to sample-representing point
    sample.color.variable=NULL,
    ellipses=T, # draw 95-CI ellipses?
    varpart.variable=NULL,
    biplot=T, # show features on beta-diversity plots?
    size.means.abundance=T, # size of feature-representing points means abundance?
    abundance.cutoff=50, # how many features to show? Deafaults to 50 most abundant ones
    feature.color.variable='Class', # which variable should be mapped to feature-representing points' color? Usually one of taxonomic levels
    triplot=F,
    chemistry.variables=NULL,
    basename="",
    analyzes=c('nmds', 'dbrda'),
    dists=c('bray', 'horn', 'GUniFrac'),
    unifrac.types=c('d_0.5'),
    statistics=T,
    find.differing.pairs=T,
    permu=999,
    p.value=0.05,
    steps=200,
    ordistep.direction='both',
    vector.col="red",
    rand.seed=667,
    type='ggplot',
    palette.family='ggsci',
    palette='nrc_npg',
    width=3.5,
    height=3.5
     ) {

	unifracs <- NULL
	disttabs <- list()

	message("plotBetaDiversity: beta-diversity analysis")
	if(!is.null(featureSet)){
	  featuretab <- featureSet$featuretab
	  sampledata <- featureSet$sampledata
	  featureannot <- featureSet$featureannot
	  tree <- featureSet$tree
	}
  set.seed(rand.seed)
  if(is.null(sample.color.variable)){
    sample.color.variable <- shape.variable
  }
  stopifnot(identical(rownames(featuretab), rownames(sampledata)), length(levels(as.factor(sampledata[[shape.variable]]))) >= 2, analyzes %in% c('nmds', 'rda', 'dbrda') );

  outlist <- list()

	for( d in dists ) {
		if( d %in% c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower","morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", "robust.aitchison")) {
			disttab <- vegan::vegdist( vegan::wisconsin(featuretab), method=d )
			disttabs[[d]] <- disttab
		} else if( grepl("UniFrac", d) ) {
			if( is.null(unifracs) ) {
				alpha <- as.numeric(sub("d_", "", unifrac.types))
				alpha <- alpha[ !is.na(alpha) ]
				unifracs <- GUniFrac::GUniFrac(featuretab, tree, alpha=alpha)$unifracs
			}
			for( a in  unifrac.types) {
				disttab <- as.dist(unifracs[ , , a ])
				name <- paste0(a, "GUniFrac")
				disttabs[[name]] <- disttab
			}
		} else {
			warning("plotBetaDiversity: unknown distance (metrics): ", d)
		}
	}

	message("plotBetaDiversity: distance matrices generated")

	outlist$disttabs <- disttabs

	outlist$analyzes <- list()

	for( analysis in analyzes ) {
			if( analysis == 'nmds' ) {
			  outlist$nmds <- list()
				for( d in names(disttabs) ){
				  message(paste0("plotBetaDiversity: ",d, " NMDS analysis"))
					outlist$analyzes$nmds[[d]] <- mgmisc1::plotNMDS(disttabs[[d]],
					                                       distname=d,
					                                       featuretab=featuretab,
					                                       sampledata=sampledata,
					                                       featureannot=featureannot,
					                                       shape.variable=shape.variable,
					                                       sample.color.variable=sample.color.variable,
					                                       feature.color.variable=feature.color.variable,
					                                       size.means.abundance=size.means.abundance,
					                                       abundance.cutoff=abundance.cutoff,
					                                       ellipses=ellipses,
					                                       biplot=biplot,
					                                       varpart.variable=varpart.variable,
  				                                       basename=basename,
					                                       statistics=statistics,
					                                       rand.seed=rand.seed,
					                                       permu=permu,
					                                       type=type,
					                                       palette.family=palette.family,
					                                       palette=palette,
					                                       width=width,
					                                       height=height)
				}
			} else if( analysis == 'dbrda' ) {
			  outlist$dbrda <- list()
				for( d in names(disttabs) ){
				  message(paste0("plotBetaDiversity: ", d, " dbRDA analysis"))
					outlist$analyzes$dbrda[[d]] <- mgmisc1::plotdbRDA(disttabs[[d]],
					                                         distname=d,
					                                         featuretab=featuretab,
					                                         sampledata=sampledata,
					                                         featureannot=featureannot,
					                                         shape.variable=shape.variable,
					                                         sample.color.variable=sample.color.variable,
					                                         varpart.variable=varpart.variable,
					                                         feature.color.variable=feature.color.variable,
					                                         size.means.abundance=size.means.abundance,
					                                         abundance.cutoff=abundance.cutoff,
					                                         ellipses=ellipses,
					                                         biplot=biplot,
					                                         triplot=triplot,
					                                         p.value=p.value,
					                                         steps=steps,
					                                         ordistep.direction=ordistep.direction,
					                                         vector.col=vector.col,
					                                         chemistry.variables=chemistry.variables,
					                                         basename=basename,
					                                         statistics=statistics,
					                                         rand.seed=rand.seed,
					                                         permu=permu,
					                                         palette.family=palette.family,
					                                         palette=palette,
					                                         type=type,
					                                         width=width,
					                                         height=height)
				}
			} else if( analysis == 'rda' ){
			  outlist$analyzes$rda <- list()
			  message("plotBetaDiversity: RDA analysis")
				outlist$analyzes$rda <- mgmisc1::plotRDA(asvtab, sampledata, shape.variable, color.variable, varpart.variable, basename=basename, statistics=statistics, type=type, palette=palette, width=width, height=height, ggplot=ggplot )
			}else if(analysis == 'cca'){
			  outlist$analyzes$cca <- list()
        message("plotBetaDiversity: CCA analysis")
        outlist$analyzes$cca <- mgmisc1::plotCCA(asvtab, sampledata, shape.variable, color.variable, basename=basename, statistics=statistics, type=type, palette=palette, width=width, height=height, ggplot=ggplot)
			}else{
			  message("plotBetaDiversity: unknown analysis ", analysis)
			}
	}
	return(outlist)
}


