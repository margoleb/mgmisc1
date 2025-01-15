#' @export palette_changer
#' @import paletteer ggsci RColorBrewer paletteer

palette_changer <- function(p, palette.family, palette, ncolors=NULL){

  if(is.null(ncolors)){
    message("palette_changer: number of colors not set")
  }else{
    message("palette_changer: number of colors: ", ncolors)
  }
  if(palette.family %in% paletteer::paletteer_packages$Name){
    palettes <- paletteer::palettes_d_names[ paletteer::palettes_d_names$package == palette.family, ]$palette
      if(palette %in% palettes){
        s <- paste0(palette.family, "::", palette)
        if(is.null(ncolors)){
          p <- p + paletteer::scale_color_paletteer_d(s)
          p <- p + paletteer::scale_fill_paletteer_d(s)
        }else{
          if(length(paletteer::paletteer_d(s)) < ncolors){
            getPalette <- colorRampPalette(paletteer::scale_color_paletteer_d(s)$palette(length(paletteer::paletteer_d(s))))
            p <- p + scale_color_manual(values=getPalette(ncolors))
            p <- p + scale_fill_manual(values=getPalette(ncolors))
          }else{
            p <- p + paletteer::scale_color_paletteer_d(s)
            p <- p + paletteer::scale_fill_paletteer_d(s)
          }
        }
      }else{
        message(paste0("palette_changer: unknown palette name for ", palette.family, " family: ", palette, "\nTo see possible palettes: paletteer::palettes_d_names"))
      }
  }else{
    message(paste0("palette_changer: unknown palette family ", palette.family, "\nTo see possible families: paletteer::paletteer_packages"))
  }

  return(p)
}
