#' @export subsetExperiment
#'

`%nin%` = Negate(`%in%`) # a negation of %in% operator, needed for samples removal to work (apparently lessR .()'s 'rows' argument can't begin with a '!')

subsetExperiment <- function(
    experiment = NULL,
    what = "samples", # either "samples" or "features"
    condition = NULL,  # a logical expression involving either variables in sampledata and their levels (for "samples"), or variables in featureannot (for "features"), must be a string (meaning it needs to be quoted)
    separate.dir = T){

  require(lessR)
  if(!is.mg(experiment)){
    stop("subsetExperiment: only an object of class 'mg' can be subsetted")
  }
  if(is.null(expr)){
    stop("subsetExperiment: a logical expression is needed")
  }

  result <- new_mg()

  result$metadata <- experiment$metadata
  subset_samples <- NULL
  subset_featureannot <- NULL
  rows <- condition
  message("rows: ", rows)

  if(what == "samples"){
    sdata <- as.data.frame(extract(experiment=experiment, feature="ASVs", what="sampledata"))
    subset_samples <- with(sdata, sdata[ eval(parse(text=rows)),])
    subset_featureannot <- extract(experiment=experiment, feature="ASVs", what="featureannot")
    seqtab <- experiment$original.data$seqtab
 #   o <- "ssamples"
 #   assign(o, subset_samples, envir = parent.frame())
#    o <- "sfeatures"
    result$original.data$sampledata <- subset_samples
    result$original.data$seqtab <- seqtab[ rownames(seqtab) %in% rownames(subset_samples), ]
    result$original.data$featureannot <- subset_featureannot
 #   assign(o, subset_features, envir = parent.frame())
  }else if(what == "features"){
    fannot <- as.data.frame(extract(experiment=experiment, feature="ASVs", what="featureannot"))
    subset_featureannot <- fannot[ .(rows), ]
    subset_samples <- extract(experiment=experiment, feature="ASVs", what="sampledata")
    subset_samples <- subset_samples[ , colnames(subset_samples) %in% rownames(subset_featureannot)]
    res$original.data <- experiment$original.data
  }else{
    stop("subsetExperiment: 'what' can be either 'samples' or 'features', and was '", what, "'")
  }

  dn <- gsub("\\ ", "_", rows)
  dn <- gsub("\\&", "and", dn)
  dn <- gsub("\\|", "or", dn)
  dn <- gsub("\\(", "", dn)
  dn <- gsub("\\)", "", dn)
  dn <- gsub("\\=\\=", "is", dn)
  dn <- gsub("\\!\\=", "is_not", dn)
  dn <- gsub("\\!is\\.null", "is_not_NULL_", dn)
  dn <- gsub("is\\.null", "is_NULL_", dn)
  dn <- gsub("\\!is\\.na", "is_not_NA_", dn)
  dn <- gsub("is\\.na", "is_NA", dn)
  dn <- gsub(">\\=", "is_greater_or_equal_to", dn)
  dn <- gsub("<\\=", "is_smaller_or_equal_to", dn)
  dn <- gsub(">", "is_greater_than", dn)
  dn <- gsub("<", "is_smaller_than", dn)


  curdir <- getwd()
  if(!is.character(separate.dir)){
    if(as.logical(separate.dir)){
      d <- paste0(curdir, "/", dn)
      message("subsetExperiment: creating directory ", d)
      if(!dir.exists(d)){
        dir.create(d)
      }
      setwd(d)
      result$metadata$working.dir <- d
    }
  }else{
    if(!dir.exists(separate.dir)){
      dir.create(separate.dir)
    }
    setwd(separate.dir)
    result$metadata$working.dir <- separate.dir
  }

  if(length(dn) > 20){
    fn <- paste0(experiment$metadata$basename, "_subset", "_mg.rds")
  }else{
    fn <- paste0(experiment$metadata$basename, "_subset", "_mg.rds")
  }

  for(f in c( "ASVs", "OTUs", "KOs", "pathways")){
    if(f == "OTUs" & !is.null(experiment$OTUs)){
      message("subsetExperiment: subsetting ", f)
      for(d in levels(experiment$OTUs$dists)){
        if(what == "features"){
          fannot <- extract(experiment=experiment, feature="OTUs", diss=d, what="featureannot")
          subset_featureannot <- fannot[ .(rows), ]
        }
        fset <- extractFeatureSet(experiment=experiment, feature="OTUs", d=d, rarefied=F)
        fset$fasta.file <- sub("\\.fasta", paste0("_", as.character(as.numeric(Sys.time())), ".fasta"), fset$fasta.file)
        fset$count.file <- sub("\\.count_table", paste0("_", as.character(as.numeric(Sys.time())), ".count_table"), fset$count.file)
        fset$tree.file <- sub("\\.tre", paste0("_", as.character(as.numeric(Sys.time())), ".tre"), fset$tree.file)
        result$OTUs$dists[[d]] <- sliceFeatureSet(set=fset, split_sampledata=subset_samples, split_featureannot=subset_featureannot)
        if(!is.null(experiment$OTUs$dists[[d]]$rarefied.data)){
          fset <- extractFeatureSet(experiment=experiment, feature="OTUs", d=d, rarefied=T)
          fset$fasta.file <- sub("\\.fasta", paste0("_", as.character(as.numeric(Sys.time())), ".fasta"), fset$fasta.file)
          fset$count.file <- sub("\\.count_table", paste0("_", as.character(as.numeric(Sys.time())), ".count_table"), fset$count.file)
          fset$tree.file <- sub("\\.tre", paste0("_", as.character(as.numeric(Sys.time())), ".tre"), fset$tree.file)
          result$OTUs$dists[[d]]$rarefied.data <- sliceFeatureSet(set=fset, split_sampledata=subset_samples, split_featureannot=subset_featureannot)
        }else{
          result$OTUs$dists[[d]]$rarefied.data <- list(NULL)
        }
      }
    }else{
      if(!is.null(experiment[[f]]$featuretab)){
        message("subsetExperiment: subsetting ", f)
        message("extracting featureset")
        fset <- extractFeatureSet(experiment=experiment, feature=f, rarefied=F)
#        o <- "fset"
#        assign(o, fset, envir=parent.frame())
        if(f == "ASVs" | what == "samples"){
          if(!is.null(fset$fasta)){
            fset$fasta.file <- sub("\\.fasta", paste0("_", as.character(as.numeric(Sys.time())), ".fasta"), fset$fasta.file)
            fset$count.file <- sub("\\.count_table", paste0("_", as.character(as.numeric(Sys.time())), ".count_table"), fset$count.file)
            fset$tree.file <- sub("\\.tre", paste0("_", as.character(as.numeric(Sys.time())), ".tre"), fset$tree.file)
          }
          message("slicing featureset")
          result[[f]] <- sliceFeatureSet(set=fset, split_sampledata=subset_samples, split_featureannot=fset$featureannot)
          message("featureset sliced")
          if(!is.null(experiment[[f]]$rarefied.data$featuretab)){
            message("extracting rarefied featureset")
            fset <- extractFeatureSet(experiment=experiment, feature=f, rarefied=T)
            fset$fasta.file <- sub("\\.fasta", paste0("_", as.character(as.numeric(Sys.time())), ".fasta"), fset$fasta.file)
            fset$count.file <- sub("\\.count_table", paste0("_", as.character(as.numeric(Sys.time())), ".count_table"), fset$count.file)
            fset$tree.file <- sub("\\.tre", paste0("_", as.character(as.numeric(Sys.time())), ".tre"), fset$tree.file)
            message("slicing rarefied featureset")
            result[[f]]$rarefied.data <- sliceFeatureSet(set=fset, split_sampledata=subset_samples, split_featureannot=fset$featureannot)
          }else{
            message("subsetExperiment: not subsetting rarefied featureset, as it is NULL")
            result[[f]]$rarefied.data <- list(NULL)
          }
        }else{
          fset <- result$ASVs

          p <- mgmisc1::preparePICRUSt2(featuretab=outlist$ASVs$featuretab,
                                        sampledata=outlist$ASVs$sampledata,
                                        fasta.file=fasta.file,
                                        filename.base=picrust.filename.base,
                                        conda=experiment$metadata$conda, condaenv=experiment$metadata$condaenv,
                                        picrust=experiment$metadata$picrust,
                                        picrust.outdir=picrust.outdir,
                                        plot.nsti=plot.nsti,
                                        variable=plot.nsti.var,
                                        processors=picrust.processors)
          if(!is.null(experiment[[f]]$rarefied.data)){
            fset <- result$ASVs$rarefied.data
          }else{
            result[[f]]$rarefied.data <- list(NULL)
          }
        }
      }else{
        message("subsetExperiment: no ", f, " to subset")
      }
    }
  }

  saveRDS(result, fn)
  return(result)
}
