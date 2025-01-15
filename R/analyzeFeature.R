#' analyzeFeature
#'
#' Analyze single feature from an mg-class experiment object.
#'
#' @param experiment An mg-class experiment object.
#' @param featurename A string indicating a name of a feature to analyze. Currently one of 'ASVs', 'OTUs', 'KOs', or 'pathways'.
#' @param d A vector of OTU dissimilarities. If featurename is 'OTUs' the argument determines which dissimilarities are analyzed. Defaults to c('0.03')
#' @param featureSet A list consisting of four elements: featuretab (matrix), sampledata (data frame), featureannot (data frame), and tree (phylo). A featureset can be extracted from an mg-class experiment with extractFeatureSet(). Disregarded if an experiment object is given.
#' @param featuretab A named integer matrix with samples in rows and features (ASVs, OTUs, genes, etc.) in columns. Can be extracted from an mg-class object with featuretab(). Disregarded if an experiment or featureset is given.
#' @param sampledata A data frame with samples in rows and environmental and grouping variables in columns. Can be extracted from an mg-class object with sampledata(). Disregarded if an experiment or featureset is given.
#' @param featureannot A data frame with features in rows and annotating variables (usually taxonomy or function) in columns.Can be extracted from an mg-class object with featureannot(). Disregarded if an experiment or featureset is given.
#' @param tree A phylo-class object or a path to a nexus file containing a tree. Can be extracted from an mg-class object with tree(). Disregarded if an mg-class object or featureset is given.
#' @param basename A string with which all filenames produced by the function will begin. Defaults to "".
#' @param directory A string being a valid path to a directory where resulting files will be stored. Defaults to ".", which means current directory.
#' @param rarefied If a featuretab/featureset is rarefied it should be set to TRUE to indicate this fact. Disregarded if an mg-class object is given. Defaults to NULL.
#' @param depth 'auto' or an integer indicating to what depth data should be rarefied. If 'auto' depth is calculated to maximize information content in rarefied data (see rarefyData()). Disregarded if an mg-class object is given or when featuretab/featureset is rarefied.
#' @param analyzes A vector of analyzes' names. Current possibilities: 'alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features'. Defaults to c('alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features').
#' @param graphics.type A string indicating if ggplot2 or base graphics should be used. Current possibilities: 'ggplot', 'base'. Defaults to 'ggplot'.
#' @param statistics A logical indicating if statistical tests should be performed. Defaults to TRUE.
#' @param permu An integer indicating a number of permutations used in statistical tests. Defaults to 999.
#' @param significance.threshold A significance threshold for tests. Defaults to 0.05.
#' @param find.differing.pairs A logical indicating if pairwise tests should be performed in case there are more than two groups and overall p-value is below significance.threshold. Defaults to TRUE.
#' @param rand.seed An odd integer being a seed for a random numbers generator. Defaults to 667 (a new better beast).
#' @param variable A string indicating name of a main variable to be analyzed. It must be one of categorical variables in sampledata. The variable is mapped to shape of symbols representing samples on beta diversity plots. Defaults to '' - must be set.
#' @param facetting.variables A vector of strings indicating maximum two variables according to which data should be divided (and displayed on plot panels). Defaults to NULL.
#' @param vertical A logical indicating if facetted plots should be organized vertically - i.e. height should be greater than width -  (TRUE) or horizontally (FALSE). Defaults to TRUE.
#' @param collect A logical indicating if common legends should be collected in one place. Defaults to FALSE.
#' @param comps A vector of strings indicating which groups should be compared. Group names must be levels of the main variable. Defaults to NULL (all possible comparisons).
#' @param processors An integer indicating number of processors to use. Defaults to 1.
#' @param palette.family A string indicating which family of palettes should be used. Any family included in the paletteer packages can be used. Defaults to 'ggsci'.
#' @param palette A string indicating a palette from palette.family. Defaults to 'nrc-npg'.
#' @param tag_level A string indicating if and how subplots of a plot should be marked. Current possibilities: NULL - no tagging, 'A' - capital letters, 'a' - small letters , '1' - arabic numbers, 'I' - capitalized Roman numbers, 'i' - small Roman numbers. Defaults to NULL.
#' @param width A number indicating width of an individual plot (in inches). Defaults to 3.5.
#' @param height A number indicating height of an individual plot (in inches). Defaults to 3.5.
#' @param alpha.div.measures A vector of strings indicating which alpha diversity measures should be used. Current possibilities: i) species richness: 'sobs' - observed richness, 'chao' - estimated total richness (Chao1 index), 'ace' - estimated total richness (ACE index); ii) diversity: 'H' - Shannon's H', 'simpson' - Simpson's D, 'invsimpson' - Simpson's 1/D; iii) evennes: 'E' - Shannon's E. Defaults to c('sobs', 'H', 'E').
#' @param alpha.plot.style A string indicating how alpha diversity should be plotted. Current possibilities: 'box', 'violin', 'dot'. Defaults to 'box'.
#' @param reflevel NULL or a string indicating a reference group for alpha diversity analysis to compare all others with. Defaults to NULL.
#' @param alpha.plot.cols An integer indicating how many columns there should be in alpha diversity plots. Defaults to 2.
#' @param alpha.subplots.annotation A string indicating if and how individual subplots of an alpha diversity plot should be tagged. Current possibilities: NULL - no tagging, 'A' - capital letters, 'a' - small letters , '1' - arabic numbers, 'I' - capitalized Roman numbers, 'i' - small Roman numbers. Defaults to NULL.
#' @param points A logical indicating if points representing individual samples should be added to alpha diversity box- or violin plots. Defaults to TRUE.
#' @param beta.diversity.analyzes A vector of strings indicating which beta diversity analyzes should be performed. Current possibilities: 'nmds', 'rda', 'dbrda', 'cca'. Defaults to c('nmds', 'dbrda')
#' @param beta.diversity.metrics A vector of strings indicating which beta diversity metrics (indices) should be used. Any metric which can be calculated with vegan's vegdist or GUniFrac. Defaults to c('bray', 'GUniFrac').
#' @param unifrac.types A vector of strings indicating which flavors of GUniFrac should be used. Current possibilities: 'd_x' where x is a float, 'd_UW' weighted UniFrac. Defaults to c('d_0.5').
#' @param sample.color.variable A string indicating which variable should be mapped to color of symbols representing samples. It must be one of categorical variables in sampledata. Defaults to ''.
#' @param sample.size.variable A string indicating which variable should be mapped to size of symbols representing samples. It must be one of continuous variables in sampledata. Defaults to ''.
#' @param ellipses A logical indicating if 95\% CI ellipses should be drawn on beta diversity plots. Defaults to TRUE.
#' @param biplot A logical indicating if biplots (plots with samples and feature) should be produced in beta diversity analyses. Defaults to TRUE.
#' @param feature.color.variable A string indicating which variable from the featureannot table should be mapped to color of symbols representing features. Defaults to ''.
#' @param size.means.abundance A logical indicating if feature abundance should be mapped to symbol size on beta diversity plots. Defaults to TRUE.
#' @param abundance.cutoff A number indicating how many features with greatest abundance should appear on a biplot. Defaults to 10.
#' @param triplot A logical indicating if a triplot should be produced in beta diversity analyses. If TRUE at least one of 'rda', 'dbrda', 'cca' needs to be included in beta.diversity.analyzes.
#' @param chemistry.variables NULL or a vector of strings indicating which variables from sampledata should be used as covariates in beta.diversity analyses. Defaults to NULL.
#' @param vector.col A string indicating color in which vectors representing significant environmental covariates should be plotted. Defaults to 'red'.
#' @param ordistep.steps An integer indicating how many steps of model selection should be performed. Defaults to 1000.
#' @param ordistep.direction A string indicating algorithm for model selection. Current possibilities: 'forward', 'backward', 'both'. Defaults to 'both'.
#' @param varpart.variable NULL or a string indicating which variable should be used in variance partitioning (to allow computing main variable's explained variance). Defaults to NULL (no variance partitioning performed).
#' @param log2FC.threshold A number indicating minimal value of abs(log2 fold change) needed to consider a feature differentially abundant. Defaults to 1 (two-fold difference).
#' @param mergeby An integer indicating by which column should results of differential analysis be joined with features annotation in featureannot. Defaults to 0 (row (i.e. feature) names).
#' @param diff.abund.methods A vector of strings with names of methods used for differential abundance analysis. Current possibilities: 'deseq' - DESeq2, 'aldex' - ALDEx2. Defaults to c('deseq').
#' @param consensus A logical indicating if only features found to be differentially abundant by all chosen methods should be regarded significant. Defaults to FALSE.
#' @param plot.diff.features A logical indicating if differentially abundant features should be plotted. Defaults to FALSE.
#' @param num.top.features An integer indicating how many significant differentially abundant features to plot. Defaults to 10.
#' @param write.single A logical indicating if tables with differentially abundant features should be written out for each comparison (pair of groups). Defaults to TRUE.
#' @param write.final A logical indicating if final table with differentially abundant features (containing all comparisons) should be written to a file. Defaults to TRUE.
#' @param mc.samples A number of MCMC samples in ALDEx2. Defaults to 128.
#' @param taxlevel.names A vector of strings with taxonomical level names. Defaults to c('phyla', 'classes', 'orders', 'families', 'genera').
#' @param rare.threshold A number between 0 and 1 indicating maximal abundance of a feature categorized as 'rare'. Defaults to 0.01
#'
#' @returns An mg-class experiment object(s) with 'analyzes' slot for a given feature populated with results of analyzes.
#'
#' @export analyzeFeature
#' @import ggplot2 ggpubr patchwork



analyzeFeature <- function(
    experiment = NULL,
    featurename = NULL,
    d = c(0.03), # OTU distance thresholds to be analyzed
    featureSet = NULL,
    featuretab = NULL,
    sampledata = NULL, # descriptions of samples (variables), samples in rows, vars in columns
    featureannot = NULL, # descriptions of features
    tree = NULL,
    basename = "",
    directory = ".",
    rarefied = NULL, # TRUE if featuretab is rarefied, otherwise FALSE, if experiment is given no need to set this var,
    depth = 'auto', # rarefaction depth, needed if a set is given, featuretable is not rarefied and feature is ASV or OTU
    analyzes = c('alpha-diversity', 'beta-diversity', 'taxonomy',  'diff-features'),
    graphics.type='ggplot',
    statistics = T,
    permu=999,
    significance.threshold=0.05,
    find.differing.pairs = T, # if overall p or q.value is significant find pairs differing significantly
    rand.seed = 667,
    variable = NULL, # main variable to be analyzed, will be used in all analyzes, needs to be categorical
    facetting.variables = NULL, # categorical variables (max. two) according to which the dataset will be divided, graphs for each unique combination of the vars will be facets of a compound plot
    vertical = T, # if the compound plot should be higher than its width (T) or wider than its height (F)
    collect = F, # if individual plots' legends should be 'collected' into one; to be
    comps = NULL, # which groups should be compared? NULL means all possible comparisons
    processors = 4,
    palette.family = 'ggsci',
    palette = 'nrc_npg',
    tag_level = NULL, # should subplots of the final plot be marked? One of 'A', 'a', '1', 'I', 'i'
    width = 3.5,
    height = 3.5,
    ### alpha-diversity analysis-specific parameters
    alpha.div.measures = c('sobs', 'H', 'E'),
    alpha.plot.style = "box", # one of "box", "violin", "dot"
    reflevel = NULL, # reference (control) group for comparisons
    alpha.plot.cols = 2,
    alpha.subplots.annotation = NULL,
    points=F, # add points to violin plots?
    ### beta-diversity analysis-specific parameters
    beta.diversity.analyzes = c('nmds', 'dbrda'),
    beta.diversity.metrics = c('bray', 'horn', 'GUniFrac'),
    unifrac.types=c('d_0.5'),
    sample.color.variable = NULL, #
    sample.size.variable = NULL, # needs to be a continuous one
    ellipses = T,
    biplot = T,
    feature.color.variable = NULL,
    size.means.abundance = T,
    abundance.cutoff = 50,
    triplot = F, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
    chemistry.variables = NULL,
    vector.col = "red",
    ordistep.steps=1000,
    ordistep.direction='both',
    varpart.variable = NULL,
    ### differentially abundant features identification-specific parameters
    log2FC.threshold = 1,
    mergeby = 0,
    diff.abund.methods = c('deseq'),
    diff.abund.taxlevels = c('phyla', 'classes', 'orders', 'families', 'genera'),
    consensus = F, # include only those features which were found to be DR by all used methods in the final results
    plot.diff.features = F, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
    num.top.features = 10,
    write.single = F,
    write.final = T, # write final, combined table to disc
    mc.samples = 128, # how many Monte Carlo samples to use in ALDEx2
    ### taxonomy analysis-specific parameters
    taxlevel.names = c('phyla', 'classes', 'orders', 'families', 'genera'),
    rare.threshold = 0.01
  ){
  result <- list()
  facets <- NULL

  all_vars <- c(variable, facetting.variables, sample.color.variable, sample.size.variable, varpart.variable, chemistry.variables)
  sdata <- sampledata(experiment=experiment, feature="ASVs", rarefied=F)
  df <- sdata[ , colnames(sdata) %in% all_vars ]
  cc <- df[ complete.cases(df), ]
  if( nrow(df) > nrow(cc) ){
    to_remove <- sdata$samplename[rownames(sdata) %nin% rownames(cc)]
    experiment <- removeSamples(experiment=experiment, samples=to_remove)
  }

  if((is.null(experiment) | !inherits(experiment, 'mg') | is.null(featurename)) & (is.null(featuretab) | is.null(sampledata)) & is.null(featureSet)){
    message("analyzeFeature: nothing to analyze, either an experiment object of class 'mg' and a feature name or a featuretab and sampledata matrices (optionally with featureannot and tree) or a featureset must be given")
    stop()
  }
  if(!is.null(experiment)){
    if(is.null(experiment[[featurename]])){
      message("analyzeFeature: no feature named '", featurename, "' in experiment object. Nothing to analyze")
      stop()
    }
  }

  if(!is.null(experiment)){
    message("analyzeFeature: extracting rarefied feature set from experiment object")
    setRarefied <- extractFeatureSet(experiment=experiment, feature=featurename, rarefied=T, d=d)
    message("analyzeFeature: extracting non-rarefied feature set from experiment object")
    setNonrarefied <- extractFeatureSet(experiment=experiment, feature=featurename, rarefied=F, d=d)
    if(is.null(experiment$working.dir)){
      experiment$working.dir <- getwd()
    }
    directory <- paste0(experiment$working.dir)
  }else if(!is.null(featuretab)){
    directory <- getwd()
    if(featurename %in% c('ASVs')){
        if(is.rarefied(featuretab)){
          setRarefied <- list( featuretab = featuretab, sampledata = sampledata, featureannot = featureannot, tree = tree)
          setNonrarefied <- NULL
          warning("analyzeFeature: rarefied featuretable, analyzing taxonomy and finding differentially represented features should be performed on a non-rarefied data")
        }else{
          setNonrarefied <-  list( featuretab = featuretab, sampledata = sampledata, featureannot = featureannot, tree = tree)
          setRarefied <- rarefyData(featuretab=featuretab, sampledata=sampledata, featureannot=featureannot, tree=tree, depth=depth)
        }
      }else{ # KOs, pathways, taxa
        if(rarefied){
          setRarefied <- list( featuretab = featuretab, sampledata = sampledata, featureannot = featureannot, tree = tree)
          setNonrarefied <- NULL
        }else{
          setNonrarefied <- list( featuretab = featuretab, sampledata = sampledata, featureannot = featureannot, tree = tree)
          setRarefied <- NULL
          warning("analyzeFeature: non-rarefied data which cannot be rarefied (", featurename, "), no alpha and beta-diversity analyzes will be done")
        }
      }
  }else{ # featureSet
    directory <- getwd()
      if(rarefied){
        setRarefied <- featureSet
        setNonrarefied <- NULL
      }else{
        setNonrarefied <- featureSet
        if(featurename %in% c('ASVs', 'OTUs')){
          setRarefied <- rarefyData(featuretab=featureSet$featuretab, sampledata=featureSet$sampledata, featureannot=featureSet$featureannot, tree=featureSet$tree, depth=depth)
        }else{
          setRarefied <- NULL
        }
      }
    }


  if(!dir.exists(directory)){
    dir.create(directory)
  }
  setwd(directory)

  for(analysis in analyzes){
    if(!dir.exists(analysis)){
      dir.create(analysis)
    }
    setwd(analysis)
    if(analysis == 'alpha-diversity'){
      if(!is.null(setRarefied)){
        if(is.null(facetting.variables)){
          bname <- paste0(basename, "_alpha.diversity_", variable)
          ### If comps=NULL all possible comparisons
          if(is.null(comps)){
            com <- combn(levels(as.factor(setRarefied$sampledata[[variable]])), 2)
            comps <- split(com, rep(1:ncol(com), each=nrow(com)))
          }
          res <- plotAlphaDiversity(    featureSet=setRarefied,
                                        featurename=featurename,
                                        variable=variable,
                                        reflevel=reflevel, # reference (control) group for comparisons
                                        annotation=alpha.subplots.annotation,
                                        measures=alpha.div.measures,
                                        basename=bname,
                                        cols=alpha.plot.cols, # number of columns in compound figures
                                        width=width,
                                        height=height,
                                        style=alpha.plot.style, # boxplot with points, dotplot, violin
                                        points=points,
                                        vjust=0.7, # distance of significance markings from brackets, the greater the parameter the smaller the distance
                                        palette.family=palette.family,
                                        palette=palette,
                                        comps=comps) # comparisons list, if NULL all possible comparisons are analyzed

            result$alpha.diversity <- res

        }else{ # there are facetting.variables, looping over their levels
          if(is.null(facets)){
            facets <- sliceData(featureSet=setRarefied, variables=facetting.variables)
            saveRDS(facets, "facets.rds")
          }
          fname <- paste0(basename, "_alpha.diversity_", variable, ".facetting_", paste(facetting.variables, collapse="_"), ".svg")
            for(f in names(facets)){
              message("analyzeFeature: analyzing alpha-diversity in ", f)
              if(length(levels(as.factor(facets[[f]]$sampledata[[variable]]))) > 1){
                message("analyzeFeature: number of levels of ", variable, ": ", length(levels(as.factor(facets[[f]]$sampledata[[variable]]))))
                bname <- paste0(basename, ".", f, "_alpha.diversity_")
                set <- list(featuretab=facets[[f]]$featuretab, sampledata=facets[[f]]$sampledata, featureannot=facets[[f]]$featureannot, tree=facets[[f]]$tree)
                facets[[f]]$alpha.diversity <- plotAlphaDiversity(
                                              featureSet=set,
                                              featurename=featurename,
                                              variable=variable,
                                              reflevel=reflevel, # reference (control) group for comparisons
                                              annotation=alpha.subplots.annotation,
                                              measures=alpha.div.measures,
                                              basename=bname,
                                              cols=alpha.plot.cols, # number of columns in compound figures
                                              width=width,
                                              height=height,
                                              style=alpha.plot.style, # boxplot with points, dotplot, violin
                                              points=points,
                                              vjust=0.7, # distance of significance markings from brackets, the greater the parameter the smaller the distance
                                              palette.family=palette.family,
                                              palette=palette,
                                              comps=comps) # comparisons list, if NULL all possible comparisons are analyzed
              }else{
                message("analyzeFeature: there is only one level of ", variable, " in ", f)
              }
            }
            p <- compose_plot(object=facets, featurename=featurename, d=d, analysis='alpha.diversity', sampledata=setRarefied$sampledata,
                              variables=facetting.variables, vertical=vertical, file.name=fname, width=width, height=height, tag_level=tag_level)
            if(is.null(result$facets)){
                result$facets <- facets
                result$alpha.diversity <- list()
                result$alpha.diversity$plot <- p
            }else{
              for(f in names(facets)){
                  result$facets[[f]]$alpha.diversity <- facets[[f]]$alpha.diversity
              }
              result$alpha.diversity <- list()
              result$alpha.diversity$plot <- p
            }
          }
        }else{
          warning("analyzeFeature: non-rarefied data and no possibility to rarefy, not performing alpha-diversity analysis")
      }
    }else if(analysis == 'beta-diversity'){
      message("analyzeFeature: analyzing beta-diversity")
      if(!is.null(setRarefied)){
        if(sum(facetting.variables %in% colnames(setRarefied$sampledata)) < length(facetting.variables)){
          stop("analyzeFeature: at least one of the facetting.variables is not a variable from sampledata\nSampledata variables: ", colnames(setRarefied$sampledata))
        }
        if(!is.null(feature.color.variable)){
          if(!feature.color.variable %in% colnames(setRarefied$featureannot)){
            warning("analyzeFeature: feature.color.variable is not a variable from the featureannot table\nNot mapping feature colors to anything\nVariables in featureannot: ", colnames(setRarefied$featureannot))
            feature.color.variable <- NULL
          }
        }
        if(is.null(facetting.variables)){
          bname <- paste0(basename)
          res <- plotBetaDiversity(
            featureSet = setRarefied,
            featurename = featurename,
            basename = bname,
            rand.seed = rand.seed,
            analyzes = beta.diversity.analyzes,
            dists = beta.diversity.metrics,
            unifrac.types = unifrac.types,
            type = graphics.type,
            palette.family = palette.family,
            palette = palette,
            width = width,
            height = height,
            statistics = statistics,
            permu = permu,
            p.value = significance.threshold,
            shape.variable = variable,
            sample.color.variable = sample.color.variable,
            ellipses = ellipses,
            varpart.variable = varpart.variable,
            biplot = biplot,
            feature.color.variable = feature.color.variable,
            size.means.abundance = size.means.abundance,
            abundance.cutoff = abundance.cutoff,
            triplot = triplot,
            chemistry.variables = chemistry.variables,
            steps = ordistep.steps,
            vector.col = vector.col
          )
          result$beta.diversity <- res
        }else{
          if(is.null(facets)){
            facets <- sliceData(featureSet=setRarefied, variables=facetting.variables)
#            f <- "facets"
#            assign(f, facets, envir = .GlobalEnv)
          }
          fname <- paste0(basename, ".facetting_", paste(facetting.variables, collapse="_")) # base name for composed graphs
          for(f in names(facets)){
            message("analyzeFeature: analyzing beta-diversity in ", f)
            if(length(levels(as.factor(facets[[f]]$sampledata[[variable]]))) > 1){
              bname <- paste0(basename, ".", f) # basename for single graphs
              set <- list(featuretab=facets[[f]]$featuretab, sampledata=facets[[f]]$sampledata, featureannot=facets[[f]]$featureannot, tree=facets[[f]]$tree)
              sr <- "setRarefied"
              assign(sr, set, envir=.GlobalEnv)
              facets[[f]]$beta.diversity <- plotBetaDiversity(
                featureSet = set,
                featurename = featurename,
                basename = bname,
                rand.seed = rand.seed,
                analyzes = beta.diversity.analyzes,
                dists = beta.diversity.metrics,
                unifrac.types = unifrac.types,
                type = graphics.type,
                palette.family = palette.family,
                palette = palette,
                width = width,
                height = height,
                statistics = statistics,
                permu = permu,
                p.value = significance.threshold,
                shape.variable = variable,
                sample.color.variable = sample.color.variable,
                ellipses = ellipses,
                varpart.variable = varpart.variable,
                biplot = biplot,
                feature.color.variable = feature.color.variable,
                size.means.abundance = size.means.abundance,
                abundance.cutoff = abundance.cutoff,
                triplot = triplot,
                chemistry.variables = chemistry.variables,
                steps = ordistep.steps,
                vector.col = vector.col
              )
            }else{
              message("analyzeFeature: there is only one level of ", variable, " in ", f)
            }
          }
#          dbg <- "facets"
#          assign(dbg, facets, envir=parent.frame())
          plots <- compose_plot(object=facets, sampledata=setRarefied$sampledata, d=d, analysis='beta.diversity', file.name=fname, variables=facetting.variables, vertical=vertical, featurename=featurename, width=width, height=height, tag_level=tag_level)

          if(!is.null(result$facets)){
            for(f in names(facets)){
              result$facets[[f]]$beta.diversity <- facets[[f]]$beta.diversity
            }
          }else{
            result$facets <- facets
          }
          result$beta.diversity <- list()
          result$beta.diversity$plots <- plots
        }
      }else{
        warning("analyzeFeature: non-rarefied data and no possibility to rarefy, not performing beta-diversity analysis")
      }
    }else if(analysis == 'taxonomy'){
      if(featurename %in% c("ASVs", "OTUs")){
        if(is.null(setNonrarefied)){
          warning("analyzeFeature: only a rarefied dataset was given, taxonomy analysis should be performed on a non-rarefied one")
          set <- setRarefied
        }else{
          set <- setNonrarefied
        }
        if(is.null(facetting.variables)){
          bname <- paste0(basename, "_taxonomy_", variable)
          res <- plotTaxonomy(
            featureSet = set,
            variable = variable,
            taxlevels = taxlevel.names,
            basename = bname,
            facetting.variables = NULL,
            rare.threshold = rare.threshold,
            vertical = vertical,
            type = graphics.type,
            palette.family = palette.family,
            palette = palette,
            height = height,
            width = width
          )
          result$taxonomy <- res
        }else{
          if(is.null(facets)){
            facets <- sliceData(featureSet=setRarefied, variables=facetting.variables)
          }
          fname <- paste0(basename, ".facetting_", paste(facetting.variables, collapse="_")) # base name for composed graphs
          for(f in names(facets)){
            message("analyzeFeature: analyzing taxonomy in ", f)
            bname <- paste0(basename, ".", f) # basename for single graphs
            set <- list(featuretab=facets[[f]]$featuretab, sampledata=facets[[f]]$sampledata, featureannot=facets[[f]]$featureannot, tree=facets[[f]]$tree)
            facets[[f]]$taxonomy <- plotTaxonomy(
              featureSet = set,
              variable = variable,
              taxlevels = taxlevel.names,
              basename = bname,
              facetting.variables = NULL,
              rare.threshold = rare.threshold,
              vertical = vertical,
              type = graphics.type,
              palette.family = palette.family,
              palette = palette,
              height = height,
              width = width
            )
          }
          if(!is.null(result$facets)){
            for(f in names(facets)){
              result$facets[[f]]$taxonomy <- facets[[f]]$taxonomy
            }
          }else{
            result$facets <- facets
          }
          # producing facetted plot
          res <- plotTaxonomy(
            featureSet = setRarefied,
            variable = variable,
            taxlevels = taxlevel.names,
            basename = fname,
            facetting.variables = facetting.variables,
            rare.threshold = rare.threshold,
            vertical = vertical,
            type = graphics.type,
            palette.family = palette.family,
            palette = palette,
            height = height,
            width = width
          )
          result$taxonomy <- res
        }
      }
    }else if( analysis == 'diff-features'){
      if(is.null(setNonrarefied)){
        warning("analyzeFeature: only a rarefied dataset was given, differential abundance analysis should be performed on a non-rarefied one")
        set <- setRarefied
      }else{
        set <- setNonrarefied
      }
      form <- as.formula(paste0("~ ", variable))
      if(is.null(facetting.variables)){
        message("analyzeFeature: finding differentially abundant features")
        bname <- paste0(basename, "_diff.features_", variable)
        res <- findDifferentiallyAbundantFeatures(
          formula = form,
          featureSet = set,
          basename = bname,
          mergeby = mergeby,
          q.threshold = significance.threshold,
          log2FC.threshold = log2FC.threshold,
          method = diff.abund.methods,
          processors = processors
        )
          result$diff.features <- list()
          result$diff.features <- res$deseq
        if(!is.null(result$taxonomy)){
          count <- 1
          for(tlevel in names(result$taxonomy$taxtabs)){
            message("analyzeFeature: finding differentially abundant taxa at the level of ", tlevel)
            bname <- paste0(basename, ".", tlevel, "_diff.features_")
            fset <- list(featuretab=result$taxonomy$taxtabs[[tlevel]], featureannot=NULL, sampledata=set$sampledata, tree=set$tree)
            res <- mgmisc1::findDifferentiallyAbundantFeatures(
              formula = form,
              featureSet = set,
              basename = bname,
              mergeby = mergeby,
              q.threshold = significance.threshold,
              log2FC.threshold = log2FC.threshold,
              method = diff.abund.methods,
              processors = processors
            )
            if(count == 1){
              finaltable <- res$deseq
            }else{
              finaltable <- rbind(finaltable, res$deseq)
            }
            result$diff.features[[tlevel]] <- res$deseq
          }
          result$diff.features$final.table <- finaltable
        }
      }else{
        if(is.null(facets)){
          facets <- mgmisc1::sliceData(featureSet=setRarefied, variables=facetting.variables)
        }
        message("analyzeFeature: finding differentially abundant features in facets")
        for(f in names(facets)){
          message("analyzeFeature: analyzing facet ", f)
          if(length(levels(as.factor(facets[[f]]$sampledata[[variable]]))) > 1){
            bname <- paste0(basename, ".", f, "_diff.features_")
            set <- list(featuretab=facets[[f]]$featuretab, sampledata=facets[[f]]$sampledata, featureannot=facets[[f]]$featureannot, tree=facets[[f]]$tree)
            facets[[f]]$diff.features <- list()
            facets[[f]]$diff.features$ASVs <- mgmisc1::findDifferentiallyAbundantFeatures(
              formula = form,
              featureSet = set,
              basename = bname,
              mergeby = mergeby,
              q.threshold = significance.threshold,
              log2FC.threshold = log2FC.threshold,
              method = diff.abund.methods,
              processors = processors
            )
            if(!is.null(result$taxonomy)){
              count <- 1
              for(tlevel in names(result$taxonomy$taxtabs)){
                if(tlevel %in% diff.abund.taxlevels){
                  message("analyzeFeature: finding differentially abundant taxa at the level of ", tlevel)
                  bname <- paste0(basename, ".", f, ".", tlevel, "_diff.features_")
                  fset <- list(featuretab=facets[[f]]$taxonomy$taxtabs[[tlevel]], featureannot=NULL, sampledata=facets[[f]]$sampledata, tree=set$tree)
                  facets[[f]]$diff.features[[tlevel]] <- mgmisc1::findDifferentiallyAbundantFeatures(
                    formula = form,
                    featureSet = fset,
                    basename = bname,
                    mergeby = mergeby,
                    q.threshold = significance.threshold,
                    log2FC.threshold = log2FC.threshold,
                    method = diff.abund.methods,
                    processors = processors
                  )$deseq
                  if(count == 1){
                    finaltable <- facets[[f]]$diff.features[[tlevel]]
                  }else{
                    finaltable <- rbind(finaltable, facets[[f]]$diff.features[[tlevel]])
                  }
                  count <- count + 1
                }
              }
              facets[[f]]$diff.features$final.table <- finaltable
            }
          }else{
            message("analyzeFeature: there is only one level of ", variable, " in ", f)
          }
        }

        if(!is.null(result$facets)){
          for(f in names(facets)){
            result$facets[[f]]$diff.features <- facets[[f]]$diff.features
          }
        }else{
          result$facets <- facets
        }
      }
    }else{
        warning("analyzeFeature: unknown analysis type: ", analysis)
    }
    setwd(directory)
  }

  return(result)
}





is.rarefied <- function(featuretab = NULL){
  res <- FALSE
  if(length(levels(as.factor(rowSums(featuretab)))) == 1){
    res <- TRUE
  }
  return(res)
}
