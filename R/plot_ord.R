#' @export plot_ord
#' @import paletteer ggplot2 ggpubr ggsci

plot_ord <- function(
  dfsamples = NULL,
  dffeatures = NULL,
  dfvectors = NULL,
  title = "",
  ellipses = T,
  xlab = "",
  ylab = "",
  shape.variable = NULL,
  sample.color.variable = NULL,
  sample.size.variable = NULL,
  feature.color.variable = NULL,
  size.means.abundance = F,
  abundance.cutoff = 50,
  angle = NULL,
  len = NULL,
  unit = NULL,
  vector.col = NULL,
  palette.family = NULL,
  palette = NULL
){
  dfsamples$shape <- as.factor(dfsamples$shape)
  dfsamples$color <- as.factor(dfsamples$color)
  if(!is.null(sample.size.variable)){
    p <- ggplot2::ggplot() + ggplot2::geom_point(data=dfsamples, aes(x=x, y=y, shape=shape.variable, color=color, size=size)) + xlab(xlab) + ylab(ylab) +
      ggplot2::theme_minimal(base_size=6) + ggplot2::theme(axis.text=element_text(size=4), legend.text=element_text(size=4)) +
      ggplot2::guides(color=guide_legend(byrow=T, title=as.character(sample.color.variable), keyheight=unit(0.2, 'cm')), shape=guide_legend(byrow=T, title=as.character(shape.variable), keyheight=unit(0.2, 'cm')), size=guide_legend(byrow=T, title=as.character(sample.size.variable), keyheight=unit(0.2, 'cm')))
  }else{
    p <- ggplot2::ggplot() + ggplot2::geom_point(data=dfsamples, aes(x=x, y=y, shape=shape, color=color)) + xlab(xlab) + ylab(ylab) +
      ggplot2::theme_minimal(base_size=6) + ggplot2::theme(axis.text=element_text(size=4), legend.text=element_text(size=4)) +
      ggplot2::guides(color=guide_legend(byrow=T, title=as.character(sample.color.variable), keyheight=unit(0.2, 'cm')), shape=guide_legend(byrow=T, title=as.character(shape.variable), keyheight=unit(0.2, 'cm')))
  }
  if(title != ""){
    p <- p + ggplot2::annotate("text", -Inf, -Inf, label=title, hjust='inward', vjust='inward', size=6/.pt) # size adjusted as per https://stackoverflow.com/questions/65076492/ggplot-size-of-annotate-vs-size-of-element-text
  }
  if(ellipses){
    p <- p + ggplot2::stat_ellipse(data=dfsamples, aes(x=x, y=y, group=shape))
  }
  if(!is.null(dffeatures)){
    if(size.means.abundance){
      p <- p + ggplot2::geom_point(data=dffeatures, aes(x=y, y=y, size=Abundance, fill=color), shape=21, color='black') +
        ggplot2::guides(fill=guide_legend(byrow=T, title=as.character(feature.color.variable), keyheight=unit(0.2, 'cm')), size=guide_legend(byrow=T, title="Abundance", keyheight=unit(0.2, 'cm')))
    }else{
      p <- p + ggplot2::geom_point(data=dffeatures, aes(x=y, y=y, fill=color), shape=21, color='black') +
        ggplot2::guides(fill=guide_legend(byrow=T, title=as.character(feature.color.variable), keyheight=unit(0.2, 'cm')))
    }
  }
  if(!is.null(dfvectors)){
    p <- p +  geom_segment(data=dfvectors, aes(x=0, xend=x, y=0, yend=y),
                           arrow=arrow(angle=angle, length=unit(len, unit)), color=vector.col) +
      geom_text(data=dfvectors, aes(x=x, y=y, label=rownames(dfvectors)), color=vector.col, hjust="outward")

    #	    p <- p + coord_fixed(ratio=1)

  }
  pal <- paste0(palette.family, "::", palette)
  p <- p + paletteer::scale_color_paletteer_d(pal)
  r <- "p"
  assign(r, p, envir=.GlobalEnv)

  return(p)
}
