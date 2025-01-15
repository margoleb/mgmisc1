#' @export plotNMDS

plotNMDS <- function(
    d, # distance matrix (usually vegdist or GUniFrac-produced)
    distname="",
    featuretab=NULL, # needed if feature (species) scores are to be computed
    sampledata=NULL,
    featureannot=NULL,
    shape.variable=NULL, # variable mapped to shape (for samples)
    sample.color.variable=NULL, # variable mapped to color (for samples)
    sample.size.variable=NULL, # variable mapped to size (for samples), will be square rooted
    varpart.variable=NULL, # additional variable used in varpart
    size.means.abundance=F, # features abundance mapped to size?
    abundance.cutoff=50, # 50 most abundant features will be ploted in biplot
    feature.color.variable=NULL, # variable mapped to color (for features)
    basename="",
    statistics=T,
    type = 'ggplot', # graphics type, either 'base' (plotting functions from the vegan package are used) or 'ggplot' (ggplot2 functions are used)
    ellipses=T, # draw 95% CI ellipses?
    biplot=F, # plot features (species)?
    palette.family='ggsci',
    palette='nrc_npg',
    permu=999,
    rand.seed=667,
    width=3.5,
    height=3.5,
    try=100,
    trymax=200){

  message("plotNMDS: plotting NMDS")
  set.seed(rand.seed)
  outlist <- list()
  outlist$varpart.variable <- varpart.variable
  outlist$shape.variable <- shape.variable
  outlist$sample.color.variable <- sample.color.variable
  outlist$sample.size.variable <- sample.size.variable
  outlist$size.means.abundance <- size.means.abundance
  outlist$feature.color.variable <- feature.color.variable
  outlist$rand.seed <- rand.seed

  explvar <- NULL # variance explained by shape variable
	permanovaf <- NULL # PERMANOVA's pseudo F
	outlist$rand.seed <- rand.seed
	outlist$varpart <- NULL
	outlist$permanova <- NULL
	ord <- vegan::metaMDS(d, k=2, try=try, trymax=trymax)
	dfsamples <- as.data.frame(vegan::scores(ord, display='sites', choices=c(1,2)))
	dfsamples$shape<- sampledata[[shape.variable]]
	dfsamples$color <- sampledata[[sample.color.variable]]
	if(!is.null(sample.size.variable)){
	  dfsamples$size <- sampledata[[sample.size.variable]]
	}
	# adding species scores
	if(biplot){
	  ord.points <- vegan::postMDS(vegan::scores(ord, display='sites', choices=c(1,2)), d)
	  dffeatures <- as.data.frame(vegan::wascores(ord.points, featuretab))
	  dffeatures$color <- featureannot[[feature.color.variable]]
	  dffeatures$Abundance <- colSums(featuretab)
	  dffeatures <- head(dffeatures[ order(dffeatures$Abundance, decreasing=T), ], abundance.cutoff)
	}
	outlist$nmds <- ord

	outlist$varpart <- NULL
	title <- ""
	if( !is.null(varpart.variable) ){
	  message("plotNMDS: performing variance partitioning")
		f1 <- as.formula(paste0(" ~ ", shape.variable))
		f2 <- as.formula(paste0(" ~ ", varpart.variable))
		ordvarpart <- vegan::varpart(d, f1, f2, data=sampledata)
		outlist$varpart <- ordvarpart
		explvar <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[1,3])), '%');
		title <- paste0("\nVariance explained by ", shape.variable,": ", explvar );
	}
	outlist$betadisper <- NULL
	outlist$betadisper.anova <- NULL
	outlist$permanova <- NULL
	if( as.logical(statistics) ) {
	  message("plotNMDS: performing PERMANOVA and betadisper analyses")
		ordbetadisper <- vegan::betadisper( d, sampledata[[shape.variable]], bias.adjust=T)
		outlist$betadisper <- ordbetadisper
		message("plotNMDS: betadisper done")
		ordbetadisperanova <- anova( ordbetadisper, permu=permu )
		outlist$betadisper.anova <- ordbetadisperanova
		message("plotNMDS: betadisper tested")
		betadisperf <- sprintf( "%.2f", ordbetadisperanova$"F value"[1] );
		betadisperp <- sprintf( "%.2e", ordbetadisperanova$"Pr(>F)"[1] );
		f <- as.formula(paste0("d ~ ", shape.variable))
		ordpermanova <- vegan::adonis2(f, data=sampledata, permu=permu)
		outlist$permanova <- ordpermanova
		message("plotNMDS: PERMANOVA done")
		permanovaf <- sprintf( "%.2f", ordpermanova$F[1] );
		permanovap <- sprintf( "%.3f", ordpermanova$"Pr(>F)"[1] );
		title <- paste0("ANOVA F = ", permanovaf, ", p = ", permanovap, "\nbetadisper F = ", betadisperf, ", p = ", betadisperp, title)
	}

	if(type == 'base'){
  	file <- paste0(basename, ".", distname, "_", shape.variable, "_", sample.color.variable, ".nmds.svg")
  	message("plotNMDS: generating plot ", file)
  	svg( file, width=width, height=height, pointsize=6 );
  	par(mar=c(4,4,4,2)+0.1, mfrow=c(1,1), las=2);
  	vegan::ordiplot(ord, type='none', main=title);
  	points( ord, pch=21+as.numeric(as.factor(sampledata[[shape.variable]])), col='black', bg=as.numeric(as.factor(sampledata[[sample.color.variable]])) );
  	vegan::ordiellipse( ord, groups=sampledata[[shape.variable]], label=T);
  	legend('topright', legend=levels(as.factor(sampledata[[shape.variable]])), title=shape.variable, pch=21+sort(unique(as.numeric(as.factor(sampledata[[shape.variable]])))));
  	legend('topleft', legend=levels(as.factor(sampledata[[sample.color.variable]])), title=sample.color.variable, pch=21, pt.bg=sort(unique(as.numeric(as.factor(sampledata[[sample.color.variable]])))));
  	dev.off();
	}else if(type =='ggplot'){
	  file <- paste0(basename, ".", distname, "_", shape.variable, "_", sample.color.variable, ".nmds.svg")
	  message("plotNMDS: generating plot ", file)
	  if(!is.null(sample.size.variable)){
	    p <- ggplot2::ggplot() + ggplot2::geom_point(data=dfsamples, aes(x=NMDS1, y=NMDS2, shape=shape, color=color, size=size)) +
	      ggplot2::theme_minimal(base_size=6) + ggplot2::theme(axis.text=element_text(size=4), legend.text=element_text(size=4)) +
	      ggplot2::guides(color=guide_legend(byrow=T, title=as.character(sample.color.variable), keyheight=unit(0.2, 'cm')), shape=guide_legend(byrow=T, title=as.character(shape.variable), keyheight=unit(0.2, 'cm')), size=guide_legend(byrow=T, title=as.character(sample.size.variable), keyheight=unit(0.2, 'cm')))

	  }else{
    p <- ggplot2::ggplot() + ggplot2::geom_point(data=dfsamples, aes(x=NMDS1, y=NMDS2, shape=shape, color=color)) +
      ggplot2::theme_minimal(base_size=6) + ggplot2::theme(axis.text=element_text(size=4), legend.text=element_text(size=4)) +
      ggplot2::guides(color=guide_legend(byrow=T, title=as.character(sample.color.variable), keyheight=unit(0.2, 'cm')), shape=guide_legend(byrow=T, title=as.character(shape.variable), keyheight=unit(0.2, 'cm')))
	  }
    if(statistics){
      p <- p + ggplot2::annotate("text", -Inf, -Inf, label=title, hjust='inward', vjust='inward', size=6/.pt) # size adjusted as per https://stackoverflow.com/questions/65076492/ggplot-size-of-annotate-vs-size-of-element-text
    }
    if(ellipses){
      p <- p + ggplot2::stat_ellipse(data=dfsamples, aes(x=NMDS1, y=NMDS2, group=shape))
    }
    if(biplot){
      if(size.means.abundance){
        p <- p + ggplot2::geom_point(data=dffeatures, aes(x=NMDS1, y=NMDS2, size=Abundance, fill=color), shape=21, color='black') +
        ggplot2::guides(fill=guide_legend(byrow=T, title=as.character(feature.color.variable), keyheight=unit(0.2, 'cm')), size=guide_legend(byrow=T, title="Abundance", keyheight=unit(0.2, 'cm')))
      }else{
        p <- p + ggplot2::geom_point(data=dffeatures, aes(x=NMDS1, y=NMDS2, fill=color), shape=21, color='black') +
          ggplot2::guides(fill=guide_legend(byrow=T, title=as.character(feature.color.variable), keyheight=unit(0.2, 'cm')))
      }
    }
    svg( file, width=width, height=height, pointsize=6 );
    print(p)
    dev.off()
    outlist$plot <- p
	}else{
	  message("plotNMDS: unknown graphics type ", type, " , no plot created")
	}
  return(outlist)
}

