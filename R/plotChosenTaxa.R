#' @export plotChosenTaxa

plotChosenTaxa <- function( asvtab, data, taxonomy=NULL, variable, reflevel=NULL, basename, taxlevel='genera', taxa=NULL, palette='default', width=3.5, height=3.5 ) {
  stopifnot( variable %in% colnames(data), identical(rownames(asvtab), rownames(data)), levels(data[[variable]]) > 1, !is.null(taxa) )
  tax <- taxonomy[ rownames(taxonomy) %in% colnames(asvtab), ]

  message("Plotting chosen taxa")

  outlist <- list()
  for( taxlevelpl in taxlevel ){
    if( taxlevelpl == "phyla" ){
      taxlevelsngl <- "Phylum";
    } else if( taxlevelpl == "classes" ){
      taxlevelsngl <- "Class";
    } else if( taxlevelpl == "orders" ){
      taxlevelsngl <- "Order";
    } else if( taxlevelpl == "families" ){
      taxlevelsngl <- "Family";
    } else {
      taxlevelsngl <- "Genus";
    }
    taxtmp <- as.data.frame(t(asvtab))
    taxtmp$taxlvl <- tax[[taxlevelsngl]];
    taxtmp <-  aggregate(. ~ taxlvl, data=taxtmp, FUN='sum');
    rownames(taxtmp) <- taxtmp$taxlvl;
    taxtmp$taxlvl <- NULL;
    taxtmp <- as.data.frame(t(taxtmp));
    outlist[[taxlevelpl]] <- taxtmp;
    message(paste0(taxlevelpl))
    taxtmp <- aggregate( . ~ data[[variable]], data=taxtmp, FUN='sum');
    rownames(taxtmp) <- paste0(taxtmp$"data[[variable]]");
    taxtmp$"data[[variable]]" <- NULL;
    taxtmp.t <- as.data.frame(t(taxtmp));
    taxtmp.t <- taxtmp.t[ order(rowSums(taxtmp.t), decreasing=T), ];
    taxtmp <- as.data.frame(t(taxtmp.t));
    taxtmp_percent <- 100 * taxtmp/rowSums(taxtmp);
    taxtmp_percent <- as.data.frame(t(taxtmp_percent))
    taxtmp_percent <- taxtmp_percent[ rownames(taxtmp_percent) %in% taxa,  ]

    if( nrow(taxtmp_percent) > 0 ){
      file <- paste0( basename, ".", taxlevelpl, ".", paste(taxa, collapse="-"),  "_", variable, ".svg" );
      palette( RColorBrewer::brewer.pal(10, "Set3") );
      svg( file, width=width, height=height, pointsize=3 );
      par( mar=c(6,4,2,12)+0.1, mfrow=c(1,1), las=2 );
      with( taxtmp_percent, barplot( as.matrix(taxtmp_percent), legend=rownames(taxtmp_percent), args.legend=list(x="right", inset=c(-0.06,0), xpd=T), col=palette()) );
      dev.off();
    } else {
      message(paste0("No chosen taxa in data"))
    }
  }
  names(outlist) <- taxlevel
  return(outlist)
}
