#' @export compose_plot
#' @import import ggplot2 patchwork

library(patchwork)

#' compose_plot
#'
#' Arrange individual plots in a grid and annotate them.
#'
#' @param object An mg-class object containing facets with 'analyzes' slot populated. Defaults to NULL
#' @param sampledata A data frame with samples as rows and environmental variables as columns. Needs to be unsplit. Defaults to NULL.
#' @param featurename A string indicating which feature is analyzed. Current possibilities: 'ASVs', 'OTUs', 'KOs', 'pathways'. Defaults to NULL.
#' @param file.name A string being either a filename or a common part of filenames (for beta diversity). Defaults to ''.
#' @param d A vector of numbers between 0 and 1 indicating which OTU dissimilarities should be processed. Needs to be given only if featurename='OTUs'. Defaults to NULL.
#' @param analysis A string indicating which type of analysis should be processed. Current possibilities: 'alpha-diversity', 'beta-diversity', 'diff-features'. Defaults to NULL.
#' @param variables A vector of strings indicating facetting variable(s). There should be one or two variables. Defaults to NULL.
#' @param vertical A logical indicating if produced plot(s) should be organized vertically (i.e. height should be greater than width). Defaults to TRUE.
#' @param collect A logical indicating if identical legends should be collected in one place. Defaults to F.
#' @param tag_level A string indicating how subplots should be tagged. Current possibilities: 'A'
#' @param width A number indicating width of a single subplot (in inches). Defaults to 3.5.
#' @param height A number indicating height of a single subplot (in inches). Defaults to 3.5.
#'
#' @returns A ggplot2 object containing composed plot or a list containing composed plots for each metric (in case of beta diversity analysis).
#' @export compose_plot
#'
#' @examples NULL
#'

compose_plot <- function(
    object = NULL,
    sampledata = NULL,
    featurename = NULL,
    file.name = "",
    d = NULL,
    analysis = NULL, # currently 'alpha-diversity', 'beta-diversity' or 'diff-features'
    variables = NULL,
    vertical = NULL,
    collect = F, # only one legend? To be used when legends of individual plots are identical
    tag_level = "A",
    width = 3.5,
    height = 3.5){

  message("compose_plot: facetting variable(s): ", variables[1], " and ", variables[2])
  p <- NULL
  v1 <- as.factor(sampledata[[variables[1]]])
  v1_length <- length(levels(as.factor(sampledata[[variables[1]]])))
  if(!is.null(variables[2])){
    v2 <- as.factor(sampledata[[variables[2]]])
    v2_length <- length(levels(as.factor(sampledata[[variables[2]]])))
  }else{
    v2_length <- 0
  }
  message("compose_plot: number of levels of first facetting variable: ", v1_length)
  if(v2_length != 0){
    message("compose_plot: number of levels of second facetting variable: ", v2_length)
  }

  n <- length(names(object))
  outlist <- list()


  if(analysis == 'beta.diversity'){
    count <- 1
    print(names(object))
    for(f in names(object)){
      print(f)
      for(a in names(object[[f]]$beta.diversity$analyses)){
        for(m in names(object[[f]]$beta.diversity$analyses[[a]])){
          message("compose_plot: facet: ", f, " analysis: ", a, " metric: ", m )
          outlist[[a]][[m]]$l <- list()
          if(!is.null(object[[f]]$beta.diversity$analyses[[a]][[m]]$plot)){
            message("compose_plot: adding plot from ", f, " analysis: ", a, " metric: ", m)
            if(count == 1){
              outlist[[a]][[m]]$p <- object[[f]]$beta.diversity$analyses[[a]][[m]]$plot
            }else if(count == 2){
              outlist[[a]][[m]]$p <- outlist[[a]][[m]]$p - object[[f]]$beta.diversity$analyses[[a]][[m]]$plot
           }else{
              outlist[[a]][[m]]$p <- outlist[[a]][[m]]$p + object[[f]]$beta.diversity$analyses[[a]][[m]]$plot
            }
          }else{
            message("compose_plot: adding an empty space, as there is no plot for ", f, " analysis: ", a, " metric: ", m)
            outlist[[a]][[m]]$p <- outlist[[a]][[m]]$p + patchwork::plot_spacer()
          }
        }
      }
      count <- count + 1
    }
  }else{
    count = 1
    p <- NULL
    l <- list()
    for(f in names(object)){
      if(!is.null(object[[f]][[analysis]]$plot)){
        message("compose_plot: adding plot from ", f)
        if(count == 1){
          p <- object[[f]][[analysis]]$plot
          l[[1]] <- p
        }else if(count == 2){
          p <- p - object[[f]][[analysis]]$plot
        }else{
          p <- p + object[[f]][[analysis]]$plot
        }
        count <- count + 1
      }else{
        message("compose_plot: adding an empty space, as there is no plot for ", f)
        p <- p + patchwork::plot_spacer()
        count <- count + 1
      }
    }
  }

  if(analysis == 'beta.diversity'){
      for(a in names(outlist)){
        for(m in names(outlist[[a]])){
          fname <- paste0(file.name, ".", a, ".", m, ".svg")
          fname1 <- paste0(file.name, ".", a, ".", m, "_wrapped.svg")
            p <- outlist[[a]][[m]]$p
            p1 <- NULL
            if(vertical){ # vertical layout
              message("compose_plot: vertical layout")
              if(v2_length != 0){
                if(v1_length >= v2_length){
                  message("compose_plot: adding plots by row")
                  if(collect){
                    p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, guides='collect', tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }else{
                    p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }
                  total.width <- v2_length * width
                  total.height <- v1_length * height
                }else{
                  message("compose_plot: adding plots by column")
                  if(collect){
                    p <- p + patchwork::plot_layout(ncol=v1_length, guides='collect', byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }else{
                    p <- p + patchwork::plot_layout(ncol=v1_length, byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }
                  total.width <- v1_length * width
                  total.height <- v2_length * height
                }
              }else{ # only one facetting variable
                if(collect){
                  p <- p + patchwork::plot_layout(nrow=v1_length, guides='collect', tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                }else{
                  p <- p + patchwork::plot_layout(nrow=v1_length, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                }
                total.width <- width
                total.height <- v1_length * height
              }
            }else{ # horizontal layout
              message("compose_plot: horizontal layout")
              if(v2_length != 0){
                if(v1_length >= v2_length){
                  message("compose_plot: adding plots by columns")
                  if(collect){
                    p <- p + patchwork::plot_layout(ncol=v1_length, guides='collect', byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }else{
                    p <- p + patchwork::plot_layout(ncol=v1_length, byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }
                  total.width <- v1_length * width
                  total.height <- v2_length * height
                }else{
                  message("compose_plot: adding plots by row")
                  if(collect){
                    p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, guides='collect', tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }else{
                    p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                  }
                  total.width <- v2_length * width
                  total.height <- v1_length * height
                }
              }else{ # only one facetting variable
                if(collect){
                  p <- p + patchwork::plot_layout(ncol=v1_length, guides='collect', tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                }else{
                  p <- p + patchwork::plot_layout(ncol=v1_length, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
                }
                total.width <- v1_length * width
                total.height <- height
              }
            }
            outlist[[a]][[m]]$p <- p
            svg(fname, width=total.width, height=total.height)
            print(p)
            dev.off()
          }
    }
  }else{
    if(vertical){ # vertical layout
      message("compose_plot: vertical layout")
      if(v2_length != 0){
        if(v1_length >= v2_length){
          message("compose_plot: adding plots by row")
          if(collect){
            p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, guides='collect', tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }else{
            p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }
          total.width <- v2_length * width
          total.height <- v1_length * height
        }else{
          message("compose_plot: adding plots by column")
          if(collect){
            p <- p + patchwork::plot_layout(ncol=v1_length, guides='collect', byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }else{
            p <- p + patchwork::plot_layout(ncol=v1_length, byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }
          total.width <- v1_length * width
          total.height <- v2_length * height
        }
      }else{ # only one facetting variable
        if(collect){
          p <- p + patchwork::plot_layout(nrow=v1_length, guides='collect', tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
        }else{
          p <- p + patchwork::plot_layout(nrow=v1_length, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
        }
        total.width <- width
        total.height <- v1_length * height
      }
    }else{ # horizontal layout
      message("compose_plot: horizontal layout")
      if(v2_length != 0){
        if(v1_length >= v2_length){
          message("compose_plot: adding plots by column")
          if(collect){
            p <- p + patchwork::plot_layout(nrow=v2_length, guides='collect', byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }else{
            p <- p + patchwork::plot_layout(nrow=v2_length, byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }
          total.width <- v1_length * width
          total.height <- v2_length * height
        }else{
          message("compose_plot: adding plots by row")
          if(collect){
            p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, guides='collect', tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }else{
            p <- p + patchwork::plot_layout(nrow=v1_length, byrow=T, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
          }
          total.width <- v2_length * width
          total.height <- v1_length * height
        }
      }else{ # only one facetting variable
        if(collect){
          p <- p + patchwork::plot_layout(ncol=v1_length, guides='collect', byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
        }else{
          p <- p + patchwork::plot_layout(ncol=v1_length, byrow=F, tag_level='new') + ggplot2::theme(legend.position='right') + patchwork::plot_annotation(tag_level=tag_level)
        }
        total.width <- v1_length * width
        total.height <- height
      }
    }
    message("compose_plot: total.width: ", total.width, ", total.height: ", total.height)
    svg(file.name, width=total.width, height=total.height)
    print(p)
    dev.off()
  }

  if(analysis == 'beta.diversity'){
    return(outlist)
  }else{
    message("compose_plot: returning single compound plot")
    return(p)
  }
}
