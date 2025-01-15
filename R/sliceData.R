#' @export sliceData

### The function returns list of slices

sliceData <- function(
    experiment = NULL,
    featureSet = NULL,
    featuretab = NULL,
    sampledata = NULL,
    tree = NULL,
    featureannot = NULL,
    variables = NULL,
    lv = NULL){ # levels to be chosen from a single variable if not NULL


  if(!is.null(featuretab) & !is.null(featureSet) & !is.null(experiment)){
    message("sliceData: either feature table, sampledata, taxonomy and tree or a featureset or an experiment object of class 'mg' may be processed")
    stop()
  }
  if(!is.null(featuretab) & is.null(sampledata)){
    message("When ASV table is processed, sampledata must be given")
    stop()
  }

  outlist <- list() # outlist is a list of 'mg' class objects

  if(is.null(experiment)){
    if(is.null(featureSet)){
      input <- list(featuretab=featuretab, sampledata=sampledata, featureannot=featureannot, tree=tree) # set for slicing
      message("sliceData: slicing data in a featuretab")
    }else{
      input <- featureSet
      message("sliceData: slicing data in a featureset")
    }
    outlist <- list() # outlist is a list containing featuretab, sampledata, taxonomy and tree
    if(is.null(lv)){
      message("sliceData: slicing hierarchically on ", paste(variables, sep=" "))
      sdata <- input$sampledata
      sdata$splitting_variables <- ""
      for(v in variables){
        sdata$splitting_variables <- paste0(sdata$splitting_variables, v, as.character(sdata[[v]]))
      }
      for(l in levels(as.factor(sdata$splitting_variables))){
        message("sliceData: slice ", l)
        split_sdata <- sdata[sdata$splitting_variables == l, ]
        outlist[[l]] <- mgmisc1::sliceFeatureSet(set=input, split_sampledata=split_sdata)
      }
    }else{
      message("sliceData: slicing featuretab samples with ", variables[1], " being ", lv, " will be chosen")
      sdata <- input$sampledata
      for(l in lv){
        split_sdata <- sdata[ sdata[[variables[1]]] == l, ]
        outlist[[l]] <- mgmisc1::sliceFeatureSet(set=input, split_sampledata=split_sdata)
      }
    }
  }else{
    experiment.name <- deparse(substitute(experiment))
    message("sliceData: slicing data in experiment object ", experiment.name)
    if(is.null(lv)){ # hierarchically splitting data according to given variables: if V1 has levels A and B and V2 has levels 1 and 2 A1, A2, B1 and B2 sets are generated
      message("sliceData: slicing hierarchically on ", paste(variables))
      message("sliceData: slicing original data in ", experiment.name)
      sdata <- mgmisc1::sampledata(experiment=experiment, feature="original.data", what="sampledata")
      sdata$splitting_variables <- ""
      for(v in variables){
        sdata$splitting_variables <- paste0(sdata$splitting_variables, v, as.character(sdata[[v]]))
      }
      experiment$original.data$sampledata$splitting_variables <- sdata$splitting_variables
      for(l in levels(as.factor(sdata$splitting_variables))){
        message("sliceData: slice ", l)
        outlist[[l]] <- mgmisc1::generateOneSlice(experiment=experiment, variable="splitting_variables", l=l)
      }
    }else{ # subset(s) containing samples with given levels of a single variable (if thera are more than one, the first is taken) are generated
      message("sliceData: slicing ", experiment, ", samples with ", variables[1], " being ", lv, " will be chosen")
      sdata <- sampledata(experiment=experiment, feature="original.data")
      for(l in lv){
        outlist[[l]] <- mgmisc1::generateOneSlice(experiment=experiment, variable=variables[1], l=l)
      }
    }
  }

  return(outlist)
}
