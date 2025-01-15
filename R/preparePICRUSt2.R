#' @export preparePICRUSt2
#' @import reticulate

preparePICRUSt2 <- function(
  featuretab = featuretab,
  sampledata = sampledata,
  variable = NULL,
  fasta.file = "seqs.fasta",
  filename.base = "picrust",
  conda = "conda", # path to conda
  condaenv = "picrust2", # name of conda environment in which PICRUSt2 is installed
  picrust = "picrust2_pipeline.py",
  picrust.outdir = "picrust2_out",
  plot.nsti = FALSE,
  processors = 1
){
  if(plot.nsti & is.null(variable)){
    stop("preparePICRUSt2: a variable must be given if an NSTI plot is to be generated")
  }
  if(dir.exists(picrust.outdir)){
    stop("preparePICRUSt2: PICRUSt2 output directory exists, either delete or specify other dir")
  }

  filename.biom = paste0(filename.base, ".biom")
  creation.date <- date()
  biom <- mgmisc1::prepareBIOM(featuretab=featuretab, filename=filename.biom, write=T)
  message("preparePICRUSt2: biom file generated")

  if(!is.null(condaenv)){
    # check if the env is already activated
    program_path <- system2("which", picrust, stdout=T)
    if(length(program_path) == 0){
      message("preparePICRUSt2: initializing conda environment in which PICRUSt2 is installed")
      reticulate::use_condaenv(conda=conda, condaenv=condaenv, required=T)
      reticulate::py_config() # python initialization, setting up env vars so that contents of conda environment's bin directory could be accessed
      message("preparePICRUSt2: ", condaenv, " conda environment activated")
    }else{
      message("preparePICRUSt2: picrust2 script not in a conda environment")
    }
  }
  picrust.version = system2(picrust, "--version")
  if( picrust.version != 0 ){
    stop("preparePICRUSt2: ", picrust, " is not a valid executable PICRUSt2 pipeline script. Maybe some conda environment should be activated?" )
  }
  message("preparePICRUSt2: now running PICRUSt2 pipeline with ", processors, " cores")
  cmd <- paste0( " -s ",  fasta.file, " -i ", filename.biom, " -o ", picrust.outdir, " -p ", processors)
  message(picrust, cmd)
  system2(picrust, cmd)

  cmd <- paste0("-i ", picrust.outdir, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO -o ", picrust.outdir, "/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz")
  system2("add_descriptions.py", cmd)

  cmd <- paste0("-i ", picrust.outdir, "/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o ", picrust.outdir, "/pathways_out/path_abun_unstrat_descrip.tsv.gz")
  system2("add_descriptions.py", cmd)

  result <- mgmisc1::readPICRUSt2Dir(dir=picrust.outdir)
  if(plot.nsti){
    filename.nstiplot <- paste0(filename.base, "_NSTI.svg")
    rownames(result$nsti) <- result$nsti$sample
    result$nsti$v<- sampledata[[variable]]
    svg(filename.nstiplot, width=7, height=7, pointsize=6)
    p <- ggplot( result$nsti, aes(x=v, y=weighted_NSTI)) + geom_boxplot() + xlab(variable) + ylab("weighted NSTI") + theme(text=element_text(size=6))
    print( p )
    dev.off()
  }
  result$sampledata <- sampledata

  return(result)
}
