#' @export plotAlphaDiversity
#' @import vegan
#' @import ggplot2
#' @import ggpubr
#' @import patchwork

plotAlphaDiversity <- function(
    featureSet=NULL,
    featurename="ASVs",
    featuretab=NULL,
    sampledata=NULL,
    variable,
    reflevel=NULL, # reference (control) group for comparisons
    annotation='a',
    measures=c('H', 'E', 'sobs'),
    basename="alpha_diversity",
    cols=2,
    width=3.5,
    height=3.5,
    style='box', # boxplot with points, dotplot, violin
    points=F,
    vjust=0.7, # distance of significance markings from brackets, the greater the parameter the smaller the distance
    palette.family='ggsci',
    palette='nrc_npg',
    comps=NULL) { # comparisons list, if NULL all possible comparisons are analyzed
  if(!is.null(featureSet)){
    featuretab <- featureSet$featuretab
    sampledata <- featureSet$sampledata
  }
	stopifnot(variable %in% colnames(sampledata), identical(rownames(featuretab), rownames(sampledata)) )
message("plotAlphaDiversity: alpha-diversity analysis")
if(length(levels(as.factor(sampledata[[variable]]))) < 2){
 message("plotAlphaDiversity: only one level of ", variable, ", nothing to plot")
  return(list())
}
	H <- NULL
	E <- NULL
	simpson <- NULL
	invsimpson <- NULL
	sobs <- NULL
	chao <- NULL
	ace <- NULL

	bdiv <- data.frame( row.names=rownames(featuretab), variable=sampledata[[variable]]);
	for( meas in measures ) {
		if( meas == 'H' ){
			bdiv$H <- vegan::diversity(featuretab)
		} else if( meas == 'simpson' ) {
			bdiv$simpson <- vegan::diversity(featuretab, index='simpson')
		} else if( meas == 'invsimpson' ) {
			bdiv$invsimpson <- vegan::diversity(featuretab, index='invsimpson')
		} else if( meas == 'E' ) {
			bdiv$E <- vegan::diversity(featuretab)/log(vegan::specnumber(featuretab))
		} else if( meas == 'sobs' ) {
			bdiv$sobs <- vegan::specnumber(featuretab)
		} else if( meas == 'chao' ) {
			bdiv$chao <- t(vegan::estimateR(featuretab))[,2]
		} else if( meas == 'ace' ) {
			bdiv$ace <- t(vegan::estimateR(featuretab))[,4]
		}
	}

  if(is.null(comps)){
    com <- combn(levels(as.factor(bdiv$variable)), 2)
    comps <- split(com, rep(1:ncol(com), each=nrow(com)))
  }

	xlabel <- gsub("_", " ", variable)

	if(nrow(featuretab) < 100){
	  point.size <- 0.4
	}else if(nrow(featuretab) < 500){
	  point.size=0.1
	}else{
	  point.size=0.01
	}


	for( meas in measures ){
		if( length(levels(as.factor(bdiv$variable))) == 2 ) {
			testtype = 'wilcox.test'
			general_p.value = NA
		} else {
			testtype = 'kruskal.test'
			general_p.value <- kruskal.test(bdiv[[meas]], bdiv$variable)$p.value
			print( paste(meas, general_p.value ) )
		}

		if( meas == 'H' ){
			title = "Diversity (Shannon's H')";
			ylab = "Shannon's H'";
			miny = min(bdiv[[meas]])
			maxy = max(bdiv[[meas]]) + 0.5 * (max(bdiv[[meas]])-min(bdiv[[meas]]))
			if(style == 'box'){
			  p_shannon <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='variable', add='jitter', add.params=list(size=point.size, color='variable', fill='variable'), width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'dot'){
			  p_shannon <- ggpubr::ggdotplot( bdiv, x = "variable", y = meas, color='variable', width=0.5, size=point.size, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'violin'){
			  p_shannon <- ggpubr::ggviolin( bdiv, x = "variable", y = meas, color='black', fill='variable', width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else{
			  stop("plotAlphaDiversity: unknown plot style. Must be one of 'box', 'dot' or 'violin'")
			}
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_shannon <- p_shannon + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=2, comparisons=comps, vjust=vjust, hide.ns=T)
			}
			p_shannon <- p_shannon + ggpubr::stat_compare_means(method=testtype, label.y.npc="bottom", label.x.npc="left", size=1.5)
		} else if( meas == 'simpson' ) {
			title = "Diversity (Simpson)"
			ylab = "Simpson's index"
			miny = floor(min(bdiv[[meas]]) - 0.1 * min(bdiv[[meas]]))
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			if(style == 'box'){
  			p_simpson <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='variable', add='jitter', add.params=list(size=point.size, color='variable', fill='variable'), width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  			  ylim(miny,maxy) + theme(legend.position="none") +
  			  ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else if(style == 'dot'){
  		  p_simpson <- ggpubr::ggdotplot( bdiv, x = "variable", y = meas, color='variable', width=0.5, size=point.size, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  		    ylim(miny,maxy) + theme(legend.position="none") +
  		    ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else if(style == 'violin'){
  		  p_simpson <- ggpubr::ggviolin( bdiv, x = "variable", y = meas, color='black', fill='variable', width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  		    ylim(miny,maxy) + theme(legend.position="none") +
  		    ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else{
  		  stop("plotAlphaDiversity: unknown plot style. Must be one of 'box', 'dot' or 'violin'")
  		}
	  if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_simpson <- p_simpson + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=2, comparisons=comps, vjust=vjust, hide.ns=T)
	  }
		p_simpson <- p_simpson  + ggpubr::stat_compare_means(method=testtype, label.y.npc="bottom", label.x.npc="left", size=1.5)
		} else if(meas == 'invsimpson') {
			title = "Diversity (inv. Simpson)"
			ylab = "Inv. Simpson's index"
			miny = floor(min(bdiv[[meas]]) - 0.1 * min(bdiv[[meas]]))
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			if(style == 'box'){
  			p_invsimpson <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='variable', add='jitter', add.params=list(size=point.size, color='variable', fill='variable'), width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  			  ylim(miny,maxy) + theme(legend.position="none") +
  			  ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else if(style == 'dot'){
  		  p_invsimpson <- ggpubr::ggdotplot( bdiv, x = "variable", y = meas, color='variable', width=0.5, size=point.size, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  		    ylim(miny,maxy) + theme(legend.position="none") +
  		    ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else if(style == 'violin'){
  		  p_invsimpson <- ggpubr::ggviolin( bdiv, x = "variable", y = meas, color='black', fill='variable', width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  		    ylim(miny,maxy) + theme(legend.position="none") +
  		    ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else{
  		  stop("plotAlphaDiversity: unknown plot style. Must be one of 'box', 'dot' or 'violin'")
  		}
	  if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_invsimpson <- p_invsimpson + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=2, comparisons=comps, vjust=vjust, hide.ns=T)
	  }
    p_invsimpson <- p_invsimpson  + ggpubr::stat_compare_means(method=testtype, label.y.npc="bottom", label.x.npc="left", size=1.5)
		} else if( meas == 'E' ){
			title = "Evenness (Shannon's E)";
			ylab = "Shannon's E";
			miny <- min(bdiv[[meas]])
			maxy <- max(bdiv[[meas]]) + 0.5 * (max(bdiv[[meas]])-min(bdiv[[meas]]))
			if(style == 'box'){
  			p_e <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='variable', add='jitter', add.params=list(size=point.size, color='variable', fill='variable'), width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  			  ylim(miny,maxy) + theme(legend.position="none") +
  			  ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else if(style == 'dot'){
  		  p_e <- ggpubr::ggdotplot( bdiv, x = "variable", y = meas, color='variable', width=0.5, size=point.size, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  		    ylim(miny,maxy) + theme(legend.position="none") +
  		    ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else if(style == 'violin'){
  		  p_e <- ggpubr::ggviolin( bdiv, x = "variable", y = meas, color='black', fill='variable', width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
  		    ylim(miny,maxy) + theme(legend.position="none") +
  		    ylab(ylab) + xlab(xlabel) + ggtitle(title)
  		}else{
  		  stop("plotAlphaDiversity: unknown plot style. Must be one of 'box', 'dot' or 'violin'")
  		}
    if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_e <- p_e + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=2, comparisons=comps, vjust=vjust, hide.ns=T)
		}
    p_e <- p_e  + ggpubr::stat_compare_means(method=testtype, label.y.npc="bottom", label.x.npc="left", size=1.5)
		} else if( meas == 'sobs'){
			title = "Species richness";
			ylab = paste0("Obs. no. ", featurename);
			miny = floor(min(bdiv[[meas]]) - 0.1 * min(bdiv[[meas]]))
			maxy = ceiling(max(bdiv[[meas]]) + 0.7 * (max(bdiv[[meas]])-min(bdiv[[meas]])))
			if(style == 'box'){
			  p_sobs <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='variable', add='jitter', add.params=list(size=point.size, color='variable', fill='variable'), width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'dot'){
			  p_sobs <- ggpubr::ggdotplot( bdiv, x = "variable", y = meas, color='variable', width=0.5, size=point.size, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'violin'){
			  p_sobs <- ggpubr::ggviolin( bdiv, x = "variable", y = meas, color='black', fill='variable', width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else{
			  stop("plotAlphaDiversity: unknown plot style. Must be one of 'box', 'dot' or 'violin'")
			}
		if( !is.na(general_p.value) & general_p.value < 0.05 ) {
			p_sobs <- p_sobs + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=2, comparisons=comps, vjust=vjust, hide.ns=T)
		}
    p_sobs <- p_sobs  + ggpubr::stat_compare_means(method=testtype, label.y.npc="bottom", label.x.npc="left", size=1.5)
		} else if( meas == 'chao' ) {
			title = "Chao1"
			ylab = paste0("Est. no. ", featurename)
			miny = floor(min(bdiv[[meas]]) - 0.1 * min(bdiv[[meas]]))
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * max(bdiv[[meas]]))
			if(style == 'box'){
			  p_chao <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='variable', add='jitter', add.params=list(size=point.size, color='variable', fill='variable'), width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'dot'){
			  p_chao <- ggpubr::ggdotplot( bdiv, x = "variable", y = meas, color='variable', width=0.5, size=point.size, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'violin'){
			  p_chao <- ggpubr::ggviolin( bdiv, x = "variable", y = meas, color='black', fill='variable', width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else{
			  stop("plotAlphaDiversity: unknown plot style. Must be one of 'box', 'dot' or 'violin'")
			}
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_chao <- p_chao + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=2, comparisons=comps, vjust=vjust, hide.ns=T)
			}
    p_chao <- p_chao  + ggpubr::stat_compare_means(method=testtype, label.y.npc="bottom", label.x.npc="left", size=1.5)
		} else if( meas == 'ace' ) {
			title = "ACE"
			ylab = paste0("Est. no. ", featurename)
			miny = floor(min(bdiv[[meas]]) - 0.1 * min(bdiv[[meas]]))
			maxy = ceiling(max(bdiv[[meas]]) + 0.5 * (max(bdiv[[meas]] - min(bdiv[[meas]]))))
			if(style == 'box'){
			  p_ace <- ggpubr::ggboxplot( bdiv, x = "variable", y = meas, color='variable', add='jitter', add.params=list(size=point.size, color='variable', fill='variable'), width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'dot'){
			  p_ace <- ggpubr::ggdotplot( bdiv, x = "variable", y = meas, color='variable', width=0.5, size=point.size, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    stat_compare_means(method=testtype, size=1.5, label.y=0.95*maxy) + ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else if(style == 'violin'){
			  p_ace <- ggpubr::ggviolin( bdiv, x = "variable", y = meas, color='black', fill='variable', width=0.5, size=0.1, ggtheme=ggpubr::theme_pubr(base_size=5), font.label=list(size=5) ) +
			    ylim(miny,maxy) + theme(legend.position="none") +
			    ylab(ylab) + xlab(xlabel) + ggtitle(title)
			}else{
			  stop("plotAlphaDiversity: unknown plot style. Must be one of 'box', 'dot' or 'violin'")
			}
			if( !is.na(general_p.value) & general_p.value < 0.05 ) {
				p_ace <- p_ace + ggpubr::stat_compare_means( method="wilcox.test", ref.group=reflevel, label="p.signif", tip.length=0.01, bracket.size=0.1, size=2, comparisons=comps, vjust=vjust, hide.ns=T)
			}
      p_ace <- p_ace  + ggpubr::stat_compare_means(method=testtype, label.y.npc="bottom", label.x.npc="left", size=1.5)
		}
	}
	f <- paste0(basename, ".", variable, ".alpha_diversity.svg")
	message("plotAlphaDiversity: generating plot ", f)
	final_plot <- (p_sobs + p_e)/(p_shannon + plot_spacer())
	if(!is.logical(annotation)){
	  final_plot <- final_plot + patchwork::plot_annotation(tag_levels=annotation)
	}else{
	  if(annotation){
	    final_plot <- final_plot + patchwork::plot_annotation(tag_levels="a")
	  }
	}
 final_plot <- final_plot + paletteer::scale_color_paletteer_d(paste0(palette.family, "::", palette))
	svg(f, width=width, height=height, pointsize=3);
	print(final_plot);
	dev.off();

  return(list(indices=bdiv, plot=final_plot))
}

