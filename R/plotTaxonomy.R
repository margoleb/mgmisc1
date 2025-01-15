#' @export plotTaxonomy
#' @import ggplot2 ggpubr ggsci ggtext reshape2

plotTaxonomy <- function(
    featureSet=NULL,
    featuretab=NULL,
    sampledata=NULL,  # sample data table for each of the rows in asvtab
    taxonomy=NULL, # taxonomy inofrmation for each of the columns in asvtab (featureannot)
    variable=NULL, # must be one of 'sampledata' columns, if not a factor, will be converted to one
    facetting.variables=NULL, # if type='ggplot' facets may be specified, first one will be in columns, second one in rows
    vertical=T, # whether facets should be vertical or horizontal
    reflevel=NULL, # reference (control) group to be on first on the left side of the graph
    rare.threshold=0.01, # minimal abundance (across all data) of a taxon not falling into the 'rare' category
    basename="bacteria", # a string, base for forming file names, taxlevel and variable names will be added after underscore, e.g. bacteria_phyla_variable.svg
    taxlevels=c('phyla', 'classes', 'orders', 'families', 'genera'), # vector of taxonomic level names in plural, if not in a classical set, final 's' will be stripped to obtain singular form
    palette.family="ggsci",
    palette='nrc_npg',
    width=3.5, # width of a single facet with up to 5 bars (length(levels(sampledata[[variable]])) <= 5)
    height=3.5, # height of a single facet
    type='ggplot', # whether use ggplot2 or base graphic (default)
    pointsize=3) {

  message("plotTaxonomy: plotting taxonomy")

if(is.null(featureSet)){
  if(is.null(featuretab) | is.null(taxonomy) | is.null(sampledata)){
    stop("plotTaxonomy: if no feature set is given, featuretab, taxonomy and sampledata need to be supplied")
  }
}else{
  featuretab <- featureSet$featuretab
  taxonomy <- featureSet$featureannot
  sampledata <- featureSet$sampledata
}

  stopifnot( variable %in% colnames(sampledata), identical(rownames(featuretab), rownames(sampledata)), length(levels(as.factor(sampledata[[variable]]))) > 1 )
  tax <- taxonomy[ rownames(taxonomy) %in% colnames(featuretab), ]

  v1_length <- length(levels(as.factor(facetting.variables[1])))
  if(!is.null(facetting.variables[2])){
    v2_length <- length(levels(as.factor(facetting.variables[2])))
  }else{
    v2_length <- 0
  }
  outlist <- list()
  plots <- list()
  base_width <- length(levels(as.factor(sampledata[[variable]]))) * (width-1)/5 + 1 # 1 inch for legend and 0.5 inch per bar if width=3.5
  base_height <- height
  message("plotTaxonomy: base facet width: ", base_width, ", base facet height: ", base_height)
  ncolors <- 0

	for( taxlevelpl in taxlevels ){
		if(taxlevelpl == 'kingdoms'){
		  taxlevelsngl <- "Kingdom"
		  } else if( taxlevelpl == "phyla" ){
			taxlevelsngl <- "Phylum";
			} else if( taxlevelpl == "classes" ){
			taxlevelsngl <- "Class";
			} else if( taxlevelpl == "orders" ){
			taxlevelsngl <- "Order";
			} else if( taxlevelpl == "families" ){
			taxlevelsngl <- "Family";
			} else if( taxlevelpl == 'genera'){
			taxlevelsngl <- "Genus";
			} else if( taxlevelpl == 'species'){
			  taxlevelsngl <- "Species"
			} else {
			  taxlevelsngl <- sub("s$", "", taxlevelpl)
			}
		taxtmp <- as.data.frame(t(featuretab))
		taxtmp$taxlvl <- taxonomy[[taxlevelsngl]];
		taxtmp <-  aggregate(. ~ taxlvl, data=taxtmp, FUN='sum');
		rownames(taxtmp) <- taxtmp$taxlvl;
		taxtmp$taxlvl <- NULL;
		taxtmp <- as.data.frame(t(taxtmp));
		outlist[[taxlevelpl]] <- taxtmp; # count matrix is returned (taxa as columns)
		message(paste0(taxlevelpl))
		taxtmp_norm <- taxtmp/rowSums(taxtmp) # normalization


    if(type == 'base'){
      vrs=c(variable)
      taxtmp <- prepare_aggregated_taxtable(taxtab=taxtmp_norm, sampledata=sampledata, vrs=vrs, rare.threshold=rare.threshold )
      rownames(taxtmp) <- paste0(taxtmp[[variable]]);
      taxtmp[[variable]] <- NULL;
      if(ncol(taxtmp) > ncolors){
        ncolors <- ncol(taxtmp)
      }
	    f <- paste0( basename, ".", taxlevelpl, "_", variable, ".svg" );
	    if(!is.null(facetting.variables)){
	      warning("plotTaxonomy: facetting variables ignored when graphics type is 'base'")
	    }
	    message("plotTaxonomy: generating plot ", f)
		  svg( f, width=base_width, height=base_height, pointsize=pointsize );
		  par( mar=c(6,4,2,12)+0.1, mfrow=c(1,1), las=2 );
		  with( taxtmp, barplot( as.matrix(t(taxtmp)), legend=colnames(taxtmp), args.legend=list(x="right", inset=c(-0.06,0), xpd=T), col=paletteer_d(paste0(palette.family, "::", palette))) );
		  dev.off();
		  p <- NULL
    } else if(type == 'ggplot'){
      vrs <- c(variable, facetting.variables)
      taxtmp_percent <- prepare_aggregated_taxtable(taxtab=taxtmp_norm, sampledata=sampledata, vrs=vrs, rare.threshold=rare.threshold)
      taxtmp_percent_molten <- reshape2::melt(taxtmp_percent, id.vars=vrs, value.name="Abundance", variable.name=taxlevelsngl)
      ntaxa <- length(levels(as.factor(taxtmp_percent_molten[[taxlevelsngl]])))
      if(ncolors < ntaxa){
        ncolors <- ntaxa
      }
      p <- ggplot(taxtmp_percent_molten, aes(x=.data[[variable]], y=Abundance, fill=.data[[taxlevelsngl]])) + geom_bar(position="fill", stat='identity') +
        ylab("Relative abundance") +
        xlab(variable) +
        labs(fill=taxlevelsngl) +
        theme_gray(base_size=8) +
        theme(axis.text=element_text(size=6), legend.text=element_text(size=6), axis.ticks.x=element_blank())
      file <- paste0( basename, ".", taxlevelpl, "_", variable, ".svg" );
      message("plotTaxonomy: filename: ", file)
      if(!is.null(facetting.variables)){
        file <- paste0(basename, ".", taxlevelpl, "_", variable, ".facets_", paste0(facetting.variables, collapse="_"), ".svg")
        if(v2_length > 0){
          message("plotTaxonomy: two facetting variables: ", facetting.variables[1], " and ", facetting.variables[2])
          if(vertical){
            if(length(levels(as.factor(sampledata[[facetting.variables[1]]]))) > length(levels(as.factor(sampledata[[facetting.variables[2]]])))){
              f <- as.formula(paste0(facetting.variables[1], " ~ ", facetting.variables[2]))
              width <- base_width * length(levels(as.factor(sampledata[[facetting.variables[2]]])))
              height <- base_height * length(levels(as.factor(sampledata[[facetting.variables[1]]])))
            }else{
              f <- as.formula(paste0(facetting.variables[2], " ~ ", facetting.variables[1]))
              width <- base_width * length(levels(as.factor(sampledata[[facetting.variables[1]]])))
              height <- base_height * length(levels(as.factor(sampledata[[facetting.variables[2]]])))
            }
          }else{
            if(length(levels(as.factor(sampledata[[facetting.variables[1]]]))) > length(levels(as.factor(sampledata[[facetting.variables[2]]])))){
              f <- as.formula(paste0(facetting.variables[2], " ~ ", facetting.variables[1]))
              width <- base_width * length(levels(as.factor(sampledata[[facetting.variables[1]]])))
              height <- base_height * length(levels(as.factor(sampledata[[facetting.variables[2]]])))
            }else{
              f <- as.formula(paste0(facetting.variables[1], " ~ ", facetting.variables[2]))
              width <- base_width * length(levels(as.factor(sampledata[[facetting.variables[2]]])))
              height <- base_height * length(levels(as.factor(sampledata[[facetting.variables[1]]])))
            }
          }
          p <- p + facet_grid(f)
        }else if(v2_length == 0){
          message("plotTaxonomy: one facetting variable: ", facetting.variables[1])
          if(vertical){
            f <- as.formula(paste0(facetting.variables, " ~ ."))
            p <- p + facet_wrap( f, nrow=length(levels(as.factor(sampledata[[facetting.variables[1]]]))) )
            width <- base_width
            height <- base_height * length(levels(as.factor(sampledata[[facetting.variables[1]]])))
            message("plotTaxonomy: width of vertical facetted plot: ", width, ", height: ", height)
          }else{
            f <- as.formula(paste0(". ~ ", facetting.variables))
            p <- p + facet_wrap( f, ncol=length(levels(as.factor(sampledata[[facetting.variables[1]]]))) )
            width <- base_width * length(levels(as.factor(sampledata[[facetting.variables[1]]])))
            height <- base_height
            message("plotTaxonomy: width of horizontal facetted plot: ", width, ", height: ", height)
          }
        }
      }

      if(palette != "default"){
        p <- palette_changer(p, palette.family, palette, ncolors)
      }
      plots[[taxlevelpl]] <- p
      message("plotTaxonomy: generating plot: ", file)
      message("plotTaxonomy: total width: ", width, ", total height: ", height)
      svg( file, width=width, height=height );
      print(p)
      dev.off()
    }
	}
  names(outlist) <- taxlevels
  result <- list(taxtabs=outlist, plots=plots)
  return(result)
}

wrap_by <- function(...) {
  ggplot2::facet_grid(ggplot2::vars(...), labeller = label_both)
}

prepare_aggregated_taxtable <- function( taxtab, sampledata, vrs, rare.threshold ){
  for(v in vrs){
    taxtab[[v]] <- sampledata[[v]]
  }
  f <- as.formula(paste(". ~ ", paste(vrs, collapse = " * ")))
  tmp <- aggregate(f, data=taxtab, FUN='mean');
  cos <- tmp[vrs]
  tmp <- tmp[ , !colnames(tmp) %in% vrs ]
  tmp.t <- as.data.frame(t(tmp));
  tmp.t <- tmp.t[ order(rowSums(tmp.t), decreasing=T), ];
  tmp <- as.data.frame(t(tmp.t));
  tmp_unclassified <- as.data.frame(tmp[ , grep("unknown|uncultured|unclassified|_ge|_fa|_or|_cl|_ph", colnames(tmp)) ]);
  tmp_classified <- as.data.frame(tmp[ , grep("unknown|uncultured|unclassified|_ge|_fa|_or|_cl|_ph", colnames(tmp), invert=T) ]);
  tmp_classified.nonrare <- as.data.frame(tmp_classified[, colSums(tmp_classified)/sum(tmp) > rare.threshold, drop=F ]);
  tmp_classified.rare <- as.data.frame(tmp_classified[, colSums(tmp_classified)/sum(tmp) <= rare.threshold]);
  tmp_classified.nonrare$rare <- rowSums(tmp_classified.rare);
  tmp_classified.nonrare$unclassifed <- rowSums(tmp_unclassified);
  tmp <- tmp_classified.nonrare;
  tmp_percent <- 100 * tmp
  tmp <- cbind(cos, as.data.frame(tmp_percent))
  return(tmp)
}
