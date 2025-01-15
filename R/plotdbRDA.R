#' @export plotdbRDA
#'
#' @import ggplot2
#' @import vegan
#' @import ggsci
#' @import ggplotify

plotdbRDA <- function(
    d, # distance matrix
    distname="", # distance name, needed for file names
    featuretab=NULL,
    sampledata=NULL,
    featureannot=NULL,
    chemistry.variables=NULL, # vector of chemistry variables' names; variables need to be in sampledata
    shape.variable,
    sample.color.variable=NULL,
    sample.size.variable=NULL,
    varpart.variable=NULL,
    ellipses=T, # plot 95 CI ellipses?
    feature.color.variable=NULL,
    biplot=F, # plot features (species)?
    size.means.abundance=T,
    abundance.cutoff=50,
    triplot=F, # plot vectors?
    p.value=0.05, # p-value threshold for continuous variables in model (those having p < p.value will be plotted)
    angle=20, # arrow angle for vectors plotting
    len=0.15, # length of arrowhead for vectors plotting
    unit='cm', # length unit
    vector.col='red', # color of plotted vectors
    basename="",
    statistics=T,
    palette.family='ggsci',
    palette='nrc_npg',
    type='base',
    sqrt.dist=F,
    rand.seed=667,
    permu=999, # number of permutations for tests
    steps=200, # max number of steps in ordistep
    ordistep.direction='both', # forward, reverse or both types of selection in ordistep?
    width=3.5,
    height=3.5) { # only for categorical variables

  set.seed(rand.seed)
  explvar <- NULL
	permanovaf <- NULL
	outlist <- list()

	dfsamples <- NULL
	dffeatures <- NULL
	dfvectors <- NULL

	message("plotdbRDA: generating dbRDA model")
	if(triplot){
	  chem <- sampledata[ , colnames(sampledata) %in% chemistry.variables ]
	  chem$shape <- sampledata[[shape.variable]]
	  # 'shotgun' model selection with ordistep
	  mod0 <- dbrda(d ~ 1, data=chem)
	  mod1 <- dbrda(d ~ ., data=chem)
	  ord <- ordistep(mod0, scope=formula(mod1), direction=ordistep.direction, permu=permu, steps=steps, Pin=p.value, Pout=2*p.value)
	  dfvectors <- as.data.frame(scores(ord, display='bp'))
	  if(length(dfvectors) > 0 ){
  	  df <- "dfvectors"
  	  assign(df, dfvectors, envir=.GlobalEnv)
  	  colnames(dfvectors) <- c('x', 'y')
  	  dfvectors <- dfvectors[ rownames(dfvectors) %in% chemistry.variables, ]
	  }else{
	    dfvectors <- NULL
	    triplot <- F
	    message("plotdbRDA: no significant vectors, falling back to unconstrained ordination")
	    f <- as.formula(paste0("d ~ ", shape.variable))
	    ord <- try(vegan::dbrda( f, data=sampledata ), silent=T)
	    if(inherits(ord, "try-error")) {
	      ord <- vegan::capscale(f, data=sampledata)
	    }
	  }
	}else{
	  f <- as.formula(paste0("d ~ ", shape.variable))
	  ord <- try(vegan::dbrda( f, data=sampledata ), silent=T)
	  if(inherits(ord, "try-error")) {
	    ord <- vegan::capscale(f, data=sampledata)
	  }
	}

	outlist$dbrda <- ord
	outlist$shape.variable <- shape.variable
	outlist$sample.color.variable <- sample.color.variable
	outlist$varpart.variable <- varpart.variable
	outlist$rand.seed <- rand.seed
	outlist$varpart <- NULL
	title <- ""
	dfsamples <- as.data.frame(vegan::scores(ord, display='sites'))
	colnames(dfsamples) <- c('x', 'y')
	dfsamples$shape <- sampledata[[shape.variable]]
	dfsamples$color <- sampledata[[sample.color.variable]]
	if(!is.null(sample.size.variable)){
	  dfsamples$size <- sampledata[[sample.size.variable]]
	  if(!is.numeric(dfsamples$size)){
	    dfsamples$size <- as.numeric(levels(as.factor(dfsamples$size)))
	  }
	}
	df <- "dfsamples"
	assign(df, dfsamples, envir = .GlobalEnv)

	if(!is.null(varpart.variable)){
	  message("plotdbRDA: Performing variance partitioning")
	  if(triplot){
	    chem.vars <- paste(chemistry.variables, collapse=" + ")
	    f1 <- as.formula(paste0(" ~ ", shape.variable))
	    f2 <- as.formula(paste0(" ~ ", chem.vars))
	    ordvarpart <- vegan::varpart(d, f1, f2, data=sampledata)
	    explvar.shape <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[1,3])), '%');
	    title <- paste0("\nVariance explained by ", shape.variable, ": ", explvar.shape );
	    explvar.chem <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[2,3])), '%')
	    outlist$varpart <- ordvarpart
	    title <- paste0(title, "\nVariance explained by physicochemical variables: ", explvar.chem );
	  }else{
  		f1 <- as.formula(paste0(" ~ ", shape.variable))
  		f2 <- as.formula(paste0(" ~ ", varpart.variable))
  		ordvarpart <- vegan::varpart(d, f1, f2, data=sampledata)
  		explvar <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[1,3])), '%');
  		title <- paste0(title, "\nVariance explained by ", shape.variable, ": ", explvar );
  		outlist$varpart <- ordvarpart
	  }
	}
	if(biplot){
	  vegan::sppscores(ord) <- sqrt(vegan::decostand(featuretab, "total"))
	  outlist$dbrda <- ord
	  dffeatures <- as.data.frame(vegan::scores(ord, display='species'))
	  colnames(dffeatures) <- c('x', 'y')
	  fan <- "feature.annot"
	  assign(fan, featureannot, envir=.GlobalEnv)
	  message("feature.color.variable: ", feature.color.variable)
    dffeatures$color <- featureannot[[feature.color.variable]]
    dffeatures$Abundance <- colSums(featuretab)
    dffeatures <- head(dffeatures[ order(dffeatures$Abundance, decreasing=T), ], abundance.cutoff)
    df <- "dffeatures"
    assign(df, dffeatures, envir=.GlobalEnv)
	}
	outlist$betadisper <- NULL
	outlist$betadisper.anova <- NULL
	outlist$permutest <- NULL
	if( as.logical(statistics) ) {
	  message("plotdbRDA: Testing model and performing betadisper analysis")
		ordbetadisper <- vegan::betadisper(as.dist(d), sampledata[[shape.variable]], bias.adjust=T) # in case of GUniFrac-produced matrices they need to be transformed to proper dist-class objects
		outlist$betadisper <- ordbetadisper
		ordbetadisperanova <- anova(ordbetadisper, permu=permu)
		outlist$betadisper.anova <- ordbetadisperanova
		betadisperf <- sprintf( "%.2f", ordbetadisperanova$"F value"[1] );
		betadisperp <- sprintf( "%.2e", ordbetadisperanova$"Pr(>F)"[1] );
		orddbrdaanova <- anova(ord, permu=permu)
		outlist$permutest <- orddbrdaanova
		dbrdaanovaf <- sprintf( "%.2f",orddbrdaanova$F[1] );
		dbrdaanovap <- sprintf( "%.3f", orddbrdaanova$"Pr(>F)"[1] );
		title <- paste0("ANOVA F = ", dbrdaanovaf, ", p = ", dbrdaanovap, "\nbetadisper F = ", betadisperf, ", p = ", betadisperp, title)
	}

	file <- paste0(basename, ".", distname, "_", shape.variable, "_", sample.color.variable, ".dbrda.svg")
	message("plotdbRDA: generating plot ", file)
	message("plotdbRDA: plot title: ", title)
	if(type == 'base'){
  	svg( file, width=width, height=height, pointsize=6 );
  	par(mar=c(4,4,4,2)+0.1, mfrow=c(1,1), las=2);
  	vegan::ordiplot(ord, type='none', main=title);
  	points( ord, pch=21+as.numeric(as.factor(sampledata[[shape.variable]])), col='black', bg=as.numeric(as.factor(sampledata[[sample.color.variable]])) );
  	vegan::ordiellipse( ord, groups=sampledata[[shape.variable]], label=T);
  	legend('topright', legend=levels(as.factor(sampledata[[shape.variable]])), title=shape.variable, pch=21+sort(unique(as.numeric(as.factor(sampledata[[shape.variable]])))));
  	legend('topleft', legend=levels(as.factor(sampledata[[sample.color.variable]])), title=sample.color.variable, pch=21, pt.bg=sort(unique(as.numeric(as.factor(sampledata[[sample.color.variable]])))));
  	dev.off();
	}else if( type == 'ggplot'){
	  message("plotdbRDA: preparing ggplot2 plot")
	  p <- plot_ord(dfsamples=dfsamples,
	                dffeatures=dffeatures,
	                dfvectors=dfvectors,
	                title=title,
	                xlab='dbRDA1',
	                ylab='dbRDA2',
	                shape.variable=shape.variable,
	                sample.color.variable=sample.color.variable,
	                feature.color.variable=feature.color.variable,
	                size.means.abundance=size.means.abundance,
	                angle=angle,
	                len=len,
	                unit=unit,
	                vector.col=vector.col,
	                palette.family=palette.family,
	                palette=palette)
	  svg( file, width=width, height=height, pointsize=6 );
	  print(p)
	  dev.off()
	  message("plotdbRDA: plot saved in ", file)
	  outlist$plot <- p
	}else{
	  message("plotdbRDA: unknown type of graphics: ", type, "; no plot generated")
	}
	return(outlist)
}


