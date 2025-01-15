#' analyzeAll
#'
#' @description
#' Analyze all features in an mg-class experiment.
#'
#' @details The function performs analyzes on various types of features resulting from one experiment. Currently the set of possible features comprises ASVs, OTUs, KO functions and pathways (from PICRUSt2 analysis). In the future other types of features will be added, such as genes or transcripts.
#'
#' @param experiment An mg-class experiment.
#' @param parameters A list or a path to a tsv file with parameters. Any other function parameter may apart from experiment may be set in this object/file. If parameters are set both in file/object, the settings can be overriden in the function call. Defaults to NULL.
#' @param dissimilarities A vector of numbers between 0 and 1 determining for which dissimilarities OTUs should be analyzed (and generated, if necessary). Defaults to c(0.03).
#' @param analyze What to analyze? A vector of strings being names of features to be analyzed. Current possibilities: 'ASVs', 'KOs', 'pathways', 'OTUs'. Defaults to c('ASVs', "OTUs', 'KOs', 'pathways').
#' @param analyzes How to analyze? A vector of strings being names of analyzes to be performed. Current possibilities: 'alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features'. Defaults to c('alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features')
#' @param processors Number of processors to use. Defaults to 1.
#' @param basename What should be a common part of produced files' names? Defaults to 'bacteria'.
#' @param outdir A path to directory where resulting plots and other files should be stored. Defaults to ".".
#' @param graphics.type Which graphics type should be used? One of 'ggplot' or 'base'. Defaults to 'ggplot'.
#' @param statistics Should statistical tests be performed? A logical. Defaults to TRUE.
#' @param rand.seed A random seed for reproducibility. Must be an odd positive integer. Defaults to 667 (a new, better beast).
#' @param significance.threshold A number between 0 and 1 determining significance threshold for statistical tests. Defaults to 0.05.
#' @param find.differing.pairs Should pairwise tests be performed to find differing group pairs if overall test is significant? A logical. Defaults to TRUE.
#' @param permu A positive integer determining the number of permutations performed during non-parametric tests such as PERMANOVA. Defaults to 999.
#' @param palette.family Which palette family should be used? Any family included in the paletteer package may be used. Defaults to 'ggsci'.
#' @param palette Which palette should be used? It should be a name of a palette from family chosen with palette.family. Any palette from the paletteer package may be used here.  Defaults to 'npg-nrc' from the 'ggsci' family.
#' @param collect A logical determining if identical legends should be collected in one place of a compound plot. Defaults to FALSE.
#' @param vertical A logical determining if compound plots should be organized vertically (i.e. their width be smaller than height). Defaults to TRUE.
#' @param width A number determining width of a single plot (in inches). Defaults to 3.5.
#' @param height A number determining height of a single plot (in inches). Defaults to 3.5.
#' @param main.variable The main variable to be analyzed. A string indicating a categorical variable from experiment's sampledata. Absolutely necessary. It will be mapped to symbol shape on beta-diversity plots and used as grouping variable in all other analyzes. Defaults to NULL.
#' @param facetting.variables Should facetted composite plots be produced? A vector of one or two variable names. Optional. The names must be chosen from colnames(sampledata) or variableNames(experiment). Defaults to NULL (no facetting).
#' @param sample.color.variable A string indicating which variable from sampledata should be mapped to color of symbols denoting samples. Optional. Must be categorical. Defaults to NULL.
#' @param sample.size.variable A string indicating which variable from sampledata should be mapped to size of symbols denoting samples. Optional. Must be continuous. Defaults to NULL.
#' @param feature.color.variable A string indicating which variable from featureannot should be mapped to color of symbols representing features. Must be categorical
#' @param chemical.variables A vector of strings indicating which variables from sampledata should be used in CCA or (db)RDA. Defaults to NULL.
#' @param varpart.variable A string indicating which variable from sampledata should be used in variance partitioning to enable assessing main variable's influence. Defaults to NULL.
#' @param alpha.div.measures Which species richness, diversity and evenness metrics should be used in alpha-diversity analysis? A vector of index names. Current possibilities: i) species richness: 'sobs' - stands for 'species observed'; 'chao' - Chao1 index (estimated total richness); 'ace' - ACE index (estimated total richness); ii) diversity: 'shannon' - Shannon's H'; 'simpson' - Simpson's D; 'invsimpson' - Simpson's 1/D; iii) evenness: 'shannoneven' - Shannon's E. Defaults to c('sobs', 'shannon', 'shannoneven').
#' @param alpha.div.points A logical determining if poits representing samples should be added to alpha diversity plots. Defaults to TRUE.
#' @param alpha.plot.style A string determining how alpha diversity will be plotted. Current possibilities: 'box', 'violin', 'dot'. Defauts to 'box'.
#' @param alpha.plot.cols A positive integer determining number of columns in compound alpha diversity plots. Defaults to 2.
#' @param beta.metrics Which beta-diversity metrics should be used? GUniFrac or any distance possible to calculate with 'vegdist' can be used. Defaults to c('bray', 'GUniFrac') - Bray-Curtis and GUniFrac.
#' @param unifracs Which variants of Generalized UniFrac should be used as beta-diversity metrics? ... Defaults to c('d_05').
#' @param beta.diversity.analyzes Which beta-diversity analyzes should be performed? Current possibilities: 'nmds', 'rda', 'dbrda'. Defaults to c('nmds', 'dbrda')
#' @param biplot A logical determining if biplots should be generated. Defaults to TRUE.
#' @param triplot A logical determining if triplots should be generated. Requires giving at least one variable in chemical.variables. Defaults to FALSE.
#' @param size.means.abundance Should size of feature symbols on a biplot convey features' abundance? Defaults to TRUE.
#' @param abundance.cutoff How many features should be displayed on a biplot? Most abundant features are displayed. An integer greater than 1. Defaults to 50.
#' @param diff.features.methods Which methods should be used to find differentially abundant features? A vector of method names. Current possibilities: 'deseq' (DESeq2 is used). Defaults to c('deseq').
#' @param log2FC.threshold A number determining a minimal absolute value of log2 fold change of a feature required to consider the feature differentially abundant. Defaults to 1 (meaning two-fold difference).
#' @param rare.threshold Threshold for identification of rare features. All features whose abundance is less than threshold will be lumped to a category named 'rare'. Must be a number between 0 and 1. Defaults to 0.01.
#'
#'
#' @returns An mg-class experiment object with 'analyzes' slots populated with results of chosen analyzes.
#'
#'
#' @export analyzeAll

analyzeAll <- function(
    # general parameters
    experiment = NULL, # an experiment object of class mg or a path to an .rds file containing one; if given, no featuretab, taxonomy, tree, design, fasta.file and count.file can be given.
    parameters = NULL,
    save = TRUE,
    analyze = c('ASVs', 'OTUs', 'KOs', 'pathways'),
    dissimilarities = c("0.03"),
    analyzes = c('alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features'),
    basename = 'bacteria',
    outdir = '.', # where results are to be stored
    graphics.type = 'ggplot', # if base ('base') graphics or ggplot2 ('ggplot') is to be used
    vertical = TRUE,
    collect = FALSE,
    tag_level = 'A',
    statistics = TRUE,
    rand.seed = 667,
    permu = 999,
    significance.threshold = 0.05,
    find.differing.pairs = TRUE,
    palette.family = 'ggsci', # palette family (package), must be one of those gathered in paletteer (https://github.com/EmilHvitfeldt/paletteer?tab=readme-ov-file), if present in the experiment object will be overwritten
    palette = 'nrc_npg', # palette name, must come from the package in 'palette.family', if present in the experiment object will be overwritten
    processors = 1,
    width = 3.5,
    height = 3.5,
    # variables
    main.variable = NULL,
    facetting.variables = NULL,
    varpart.variable = NULL,
    sample.color.variable = NULL,
    sample.size.variable = NULL,
    feature.color.variable = NULL,
    chemistry.variables = NULL,
    # alpha diversity
    alpha.div.measures = c('sobs', 'H', 'E'),
    alpha.div.points = T,
    ref.group = NULL,
    alpha.plot.style = 'box',
    alpha.plot.cols = 2,
    # beta diversity
    beta.diversity.metrics = c('bray', 'GUniFrac'), # one or more of metrics implemented in vegan::vegdist or GUniFrac
    unifrac.types = c('d_0.5'), # GUniFrac metrics (for ASVs and OTUs only)
    beta.diversity.analyzes = c('nmds', 'dbrda'),
    biplot=T,
    triplot=F,
    size.means.abundance=T,
    abundance.cutoff=50,
    ellipses = T,
    vector.col = 'red',
    ordistep.steps=500,
    ordistep.direction='both',
    # differentially abundant features identification parameters
    diff.features.methods=c('deseq'), # one of 'deseq', 'aldex',
    log2FC.threshold = 1,
    mergeby = 0,
    consensus = F, # include only those features which were found to be DR by all used methods in the final results
    plot.diff.features = F, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
    num.top.features = 10,
    write.single = T,
    write.final = T, # write final, combined table to disc
    mc.samples = 128, # how many Monte Carlo samples to use in ALDEx2
    # taxonomy plots parameters
    rare.threshold = 0.01,
    taxlevel.names = c('kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera')
){
  argv <- as.list(match.call()[-1])
  params <- set_parameters(p=parameters, argv)

  set.seed(params$rand.seed)

  if(is.null(experiment)){
    stop("analyzeAll: an object of class 'mg' needs to be given")
  }

  if(!R.utils::isAbsolutePath(params$outdir)){
    outdir <- file.path(getwd(), params$outdir)
  }
  if(!dir.exists(outdir)){
    dir.create(outdir)
  }
  setwd(outdir)

  result <- experiment

  if(!is.null(experiment)){
    if(is.character(experiment)){
      experiment <- readRDS(experiment)
    }
    if(!inherits(experiment, "mg")){
      stop("analyzeAll: an experiment to analyze must be of class 'mg'")
    }
    if(!is.null(params$palette.family)){
      experiment$metadata$palette.family <- params$palette.family
    }else if(!is.null(experiment$metadata$palette.family)){
      palette.family <- experiment$metadata$palette.family
    }else{
      palette.family <- 'ggsci'
    }
    if(!is.null(params$palette)){
      experiment$metadata$palette <- params$palette
    }else if(!is.null(experiment$metadata$palette)){
      palette <- experiment$metadata$palette
    }else{
      palette <- 'nrc_npg'
    }

    dissimilarities <- as.character(params$dissimilarities)
    # check if OTUs and PICRUSt data are generated, generate if needed
    experiment <- check_and_generate_data(experiment=experiment, analyze=params$analyze, dissimilarities=dissimilarities, analyzes=analyzes, processors=processors)

    for( feature in params$analyze ){
      if( feature == 'ASVs'){
        message("analyzeAll: analyzing ASVs")
        d <- paste0("ASVs_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
        if(!dir.exists(d)){
          dir.create(d)
        }
        bn <- paste0(params$basename, "_ASVs")
        message("analyzeAll: results will be stored in the ", d)
        setwd(d)
        result$ASVs$analyzes <- analyzeFeature(
          experiment=experiment,
          featurename="ASVs",
          basename=bn,
          directory=".",
          analyzes=params$analyzes,
          graphics.type=params$graphics.type,
          vertical=params$vertical,
          collect=params$collect,
          palette.family=params$palette.family,
          palette=params$palette,
          tag_level=params$tag_level,
          width=params$width,
          height=params$height,
          statistics=params$statistics,
          permu=params$permu,
          significance.threshold=params$significance.threshold,
          find.differing.pairs=params$find.differing.pairs,
          rand.seed=params$rand.seed,
          variable=params$main.variable,
          facetting.variables=params$facetting.variables,
          alpha.div.measures=params$alpha.div.measures,
          alpha.plot.style=params$alpha.plot.style,
          alpha.plot.cols=params$alpha.plot.cols,
          reflevel=params$ref.group,
          alpha.subplots.annotation=params$alpha.subplots.annotation,
          points=params$alpha.div.points,
          beta.diversity.analyzes=params$beta.diversity.analyzes,
          beta.diversity.metrics=params$beta.diversity.metrics,
          unifrac.types=params$unifrac.types,
          sample.color.variable = params$sample.color.variables, #
          sample.size.variable = params$sample.size.variable, # needs to be a continuous one
          ellipses = params$ellipses,
          biplot = params$biplot,
          feature.color.variable = params$feature.color.variable,
          size.means.abundance = params$size.means.abundance,
          abundance.cutoff = params$abundance.cutoff,
          triplot = params$triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
          chemistry.variables = params$chemical.variables,
          vector.col = params$vector.col,
          ordistep.steps=params$ordistep.steps,
          ordistep.direction=params$ordistep.direction,
          varpart.variable = params$varpart.variable,
          ### differentially abundant features identification-specific parameters
          log2FC.threshold = params$log2FC.threshold,
          mergeby = params$mergeby,
          diff.abund.methods = params$diff.features.methods,
          consensus = params$consensus, # include only those features which were found to be DR by all used methods in the final results
          plot.diff.features = params$plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
          num.top.features = params$num.top.features,
          write.single = params$write.single,
          write.final = params$write.final, # write final, combined table to disc
          mc.samples = params$mc.samples, # how many Monte Carlo samples to use in ALDEx2
          ### taxonomy analysis-specific parameters
          taxlevel.names = params$taxlevel.names,
          rare.threshold = params$rare.threshold
        )
        setwd(outdir)
        if(save){
          name <- paste0(experiment$metadata$basename, "_mg.rds")
          saveRDS(result, name)
        }
      }else if( feature == 'OTUs' ){
        message("analyzeAll: analyzing OTUs")
        for(diss in names(experiment$OTUs$dists)){
          message("analyzeAll: analyzing OTUs at ", diss, " dissimilarity")
          d <- paste0("OTUs_", diss, "_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
          dir.create(d)
          message("analyzeAll: results will be stored in the ", d)
          setwd(d)
          bn <- paste0(basename, "_OTUs")
          result$OTUs$dists[[diss]]$analyzes <- analyzeFeature(
            experiment=experiment,
            featurename="OTUs",
            d=diss,
            basename=bn,
            directory=".",
            analyzes=params$analyzes,
            graphics.type=params$graphics.type,
            vertical=params$vertical,
            collect=params$collect,
            palette.family=params$palette.family,
            palette=params$palette,
            tag_level=params$tag_level,
            width=params$width,
            height=params$height,
            statistics=params$statistics,
            permu=params$permu,
            significance.threshold=params$significance.threshold,
            find.differing.pairs=params$find.differing.pairs,
            rand.seed=params$rand.seed,
            variable=params$main.variable,
            facetting.variables=params$facetting.variables,
            alpha.div.measures=params$alpha.div.measures,
            alpha.plot.style=params$alpha.plot.style,
            alpha.plot.cols=params$alpha.plot.cols,
            reflevel=params$ref.group,
            alpha.subplots.annotation=params$alpha.subplots.annotation,
            points=params$alpha.div.points,
            beta.diversity.analyzes=params$beta.diversity.analyzes,
            beta.diversity.metrics=params$beta.diversity.metrics,
            unifrac.types=params$unifrac.types,
            sample.color.variable = params$sample.color.variables, #
            sample.size.variable = params$sample.size.variable, # needs to be a continuous one
            ellipses = params$ellipses,
            biplot = params$biplot,
            feature.color.variable = params$feature.color.variable,
            size.means.abundance = params$size.means.abundance,
            abundance.cutoff = params$abundance.cutoff,
            triplot = params$triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
            chemistry.variables = params$chemical.variables,
            vector.col = params$vector.col,
            ordistep.steps=params$ordistep.steps,
            ordistep.direction=params$ordistep.direction,
            varpart.variable = params$varpart.variable,
            ### differentially abundant features identification-specific parameters
            log2FC.threshold = params$log2FC.threshold,
            mergeby = params$mergeby,
            diff.abund.methods = params$diff.features.methods,
            consensus = params$consensus, # include only those features which were found to be DR by all used methods in the final results
            plot.diff.features = params$plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
            num.top.features = params$num.top.features,
            write.single = params$write.single,
            write.final = params$write.final, # write final, combined table to disc
            mc.samples = params$mc.samples, # how many Monte Carlo samples to use in ALDEx2
            ### taxonomy analysis-specific parameters
            taxlevel.names = params$taxlevel.names,
            rare.threshold = params$rare.threshold
          )
          setwd(outdir)
          if(save){
            name <- paste0(experiment$metadata$basename, "_mg.rds")
            saveRDS(result, name)
          }
        }
      }else if( feature == 'KOs'){
        message("analyzeAll: analyzing KOs")
        d <- paste0("KOs_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
        dir.create(d)
        message("analyzeAll: results will be stored in the ", d)
        setwd(d)
        bn <- paste0(basename, "_KOs")
        bdm <- beta.diversity.metrics[ -grep("GUniFrac", beta.diversity.metrics) ]
        result$KOs$analyzes <- analyzeFeature(
          experiment=experiment,
          featurename="KOs",
          d=NULL,
          basename=bn,
          directory=".",
          analyzes=params$analyzes,
          graphics.type=params$graphics.type,
          vertical=params$vertical,
          collect=params$collect,
          palette.family=params$palette.family,
          palette=params$palette,
          tag_level=params$tag_level,
          width=params$width,
          height=params$height,
          statistics=params$statistics,
          permu=params$permu,
          significance.threshold=params$significance.threshold,
          find.differing.pairs=params$find.differing.pairs,
          rand.seed=params$rand.seed,
          variable=params$main.variable,
          facetting.variables=params$facetting.variables,
          alpha.div.measures=params$alpha.div.measures,
          alpha.plot.style=params$alpha.plot.style,
          alpha.plot.cols=params$alpha.plot.cols,
          reflevel=params$ref.group,
          alpha.subplots.annotation=params$alpha.subplots.annotation,
          points=params$alpha.div.points,
          beta.diversity.analyzes=params$beta.diversity.analyzes,
          beta.diversity.metrics=bdm,
          unifrac.types=NULL,
          sample.color.variable = params$sample.color.variables, #
          sample.size.variable = params$sample.size.variable, # needs to be a continuous one
          ellipses = params$ellipses,
          biplot = F,
          feature.color.variable = NULL,
          size.means.abundance = params$size.means.abundance,
          abundance.cutoff = params$abundance.cutoff,
          triplot = params$triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
          chemistry.variables = params$chemical.variables,
          vector.col = params$vector.col,
          ordistep.steps=params$ordistep.steps,
          ordistep.direction=params$ordistep.direction,
          varpart.variable = params$varpart.variable,
          ### differentially abundant features identification-specific parameters
          log2FC.threshold = params$log2FC.threshold,
          mergeby = params$mergeby,
          diff.abund.methods = params$diff.features.methods,
          diff.abund.taxlevels = params$diff.abund.taxlevels,
          consensus = params$consensus, # include only those features which were found to be DR by all used methods in the final results
          plot.diff.features = params$plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
          num.top.features = params$num.top.features,
          write.single = params$write.single,
          write.final = params$write.final, # write final, combined table to disc
          mc.samples = params$mc.samples, # how many Monte Carlo samples to use in ALDEx2
          ### taxonomy analysis-specific parameters
          taxlevel.names = params$taxlevel.names,
          rare.threshold = params$rare.threshold
        )
        setwd(outdir)
        if(save){
          name <- paste0(experiment$metadata$basename, "_mg.rds")
          saveRDS(result, name)
        }
      }else if( feature == 'pathways'){
        message("analyzeAll: analyzing pathways")
        d <- paste0("pathways_", gsub(":", "-", gsub("\\ ", "_", as.character(date()))))
        dir.create(d)
        message("analyzeAll: results will be stored in ", d)
        setwd(d)
        bn <- paste0(basename, "_pathways")
        bdm <- beta.diversity.metrics[ -grep("GUniFrac", beta.diversity.metrics) ]
        a <- params$analyzes[-grep("taxonomy", params$analyzes)]
        result$pathways$analyzes <- analyzeFeature(
          experiment=experiment,
          featurename="pathways",
          d=NULL,
          basename=bn,
          directory=".",
          analyzes=a,
          graphics.type=params$graphics.type,
          vertical=params$vertical,
          collect=params$collect,
          palette.family=params$palette.family,
          palette=params$palette,
          tag_level=params$tag_level,
          width=params$width,
          height=params$height,
          statistics=params$statistics,
          permu=params$permu,
          significance.threshold=params$significance.threshold,
          find.differing.pairs=params$find.differing.pairs,
          rand.seed=params$rand.seed,
          variable=params$main.variable,
          facetting.variables=params$facetting.variables,
          alpha.div.measures=params$alpha.div.measures,
          alpha.plot.style=params$alpha.plot.style,
          alpha.plot.cols=params$alpha.plot.cols,
          reflevel=params$ref.group,
          alpha.subplots.annotation=params$alpha.subplots.annotation,
          points=params$alpha.div.points,
          beta.diversity.analyzes=params$beta.diversity.analyzes,
          beta.diversity.metrics=bdm,
          unifrac.types=NULL,
          sample.color.variable = params$sample.color.variables, #
          sample.size.variable = params$sample.size.variable, # needs to be a continuous one
          ellipses = params$ellipses,
          biplot = F,
          feature.color.variable = NULL,
          size.means.abundance = params$size.means.abundance,
          abundance.cutoff = params$abundance.cutoff,
          triplot = params$triplot, #' plot vectors representing chemical variables significantly influencing community structure? For this to take effect 'rda', 'dbrda' or 'cca' must be in beta.diversity.analyzes
          chemistry.variables = params$chemical.variables,
          vector.col = params$vector.col,
          ordistep.steps=params$ordistep.steps,
          ordistep.direction=params$ordistep.direction,
          varpart.variable = params$varpart.variable,
          ### differentially abundant features identification-specific parameters
          log2FC.threshold = params$log2FC.threshold,
          mergeby = params$mergeby,
          diff.abund.methods = params$diff.features.methods,
          consensus = params$consensus, # include only those features which were found to be DR by all used methods in the final results
          plot.diff.features = params$plot.diff.features, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
          num.top.features = params$num.top.features,
          write.single = params$write.single,
          write.final = params$write.final, # write final, combined table to disc
          mc.samples = params$mc.samples, # how many Monte Carlo samples to use in ALDEx2
          ### taxonomy analysis-specific parameters
          taxlevel.names = params$taxlevel.names,
          rare.threshold = params$rare.threshold
        )
        setwd(outdir)
        if(save){
          name <- paste0(experiment$metadata$basename, "_mg.rds")
          saveRDS(result, name)
        }
      }
      setwd(outdir)
    }
  }

  return(result)
}


set_parameters <- function(
    p=parameters, # a path to a tsv file with two columns: par_names and values or a named list where names are parameter names and elements are values
    a=argv){  # list of arguments of a calling function

  # The function returns list of parameters for analyzeAll, if a param was defined in a parameters file or in a command line it is set to the modified value,
  # otherwise default is returned. Command line has precedence over params file (i.e. if a param was defined in a file and in a command line, command line is
  # returned).

  defaults <- list(
    save = TRUE,
    analyze = c('ASVs', 'OTUs', 'KOs', 'pathways'),
    dissimilarities = c("0.03"),
    analyzes = c('alpha-diversity', 'beta-diversity', 'taxonomy', 'diff-features'),
    basename = 'bacteria',
    outdir = '.', # where results are to be stored
    graphics.type = 'ggplot', # if base ('base') graphics or ggplot2 ('ggplot') is to be used
    vertical = TRUE,
    collect = FALSE,
    tag_level = 'A',
    statistics = TRUE,
    rand.seed = 667,
    permu = 999,
    significance.threshold = 0.05,
    find.differing.pairs = TRUE,
    palette.family = 'ggsci', # palette family (package), must be one of those gathered in paletteer (https://github.com/EmilHvitfeldt/paletteer?tab=readme-ov-file), if present in the experiment object will be overwritten
    palette = 'nrc_npg', # palette name, must come from the package in 'palette.family', if present in the experiment object will be overwritten
    processors = 1,
    width = 3.5,
    height = 3.5,
    # variables
    main.variable = NULL,
    facetting.variables = NULL,
    varpart.variable = NULL,
    sample.color.variable = NULL,
    sample.size.variable = NULL,
    feature.color.variable = NULL,
    chemistry.variables = NULL,
    # alpha diversity
    alpha.div.measures = c('sobs', 'H', 'E'),
    alpha.div.points = T,
    ref.group = NULL,
    alpha.plot.style = 'box',
    alpha.plot.cols = 2,
    # beta diversity
    beta.diversity.metrics = c('bray', 'GUniFrac'), # one or more of metrics implemented in vegan::vegdist or GUniFrac
    unifrac.types = c('d_0.5'), # GUniFrac metrics (for ASVs and OTUs only)
    beta.diversity.analyzes = c('nmds', 'dbrda'),
    biplot=T,
    triplot=F,
    size.means.abundance=T,
    abundance.cutoff=50,
    ellipses = T,
    vector.col = 'red',
    ordistep.steps=500,
    ordistep.direction='both',
    # differentially abundant features identification parameters
    diff.features.methods=c('deseq'), # one of 'deseq', 'aldex',
    diff.abund.taxlevels=c('phyla', 'classes', 'orders', 'families', 'genera'),
    log2FC.threshold = 1,
    mergeby = 0,
    consensus = F, # include only those features which were found to be DR by all used methods in the final results
    plot.diff.features = F, # makes sense if the number of groups in variable is 2 or three (then a ternary plot is produced), otherwise many plots with comps will be produced
    num.top.features = 10,
    write.single = T,
    write.final = T, # write final, combined table to disc
    mc.samples = 128, # how many Monte Carlo samples to use in ALDEx2
    # taxonomy plots parameters
    rare.threshold = 0.01,
    taxlevel.names = c('kingdoms', 'phyla', 'classes', 'orders', 'families', 'genera')
  )


  if(is.character(p)){
    pdf <- read.table(p, header=T, sep="\t", dec=".", strip.white=T)
    if(!identical(colnames(pdf), c("par_names", "values"))){
      stop(paste0("set_parameters: Incorrectly formatted paramters file ", p, ". Columns should be named 'par_names' and 'values'\n"))
    }
    p <- list()
    for(n in 1:nrow(pdf)){
      p[[n]] <- pdf[n,2]
    }
    names(p) <- pdf[,1]
  }
  for(n in names(p)){
    p[[n]] <- strsplit(p[[n]], ",")[[1]]
    if(n %in% names(defaults)){
        print(paste(n, defaults[[n]], p[[n]]))
        if(is.logical(defaults[[n]])){
          defaults[[n]] <- as.logical(p[[n]])
        }else if(is.numeric(defaults[[n]])){
          defaults[[n]] <- as.numeric(p[[n]])
        }else{
          defaults[[n]] <- p[[n]]
        }
      if(length(defaults[[n]]) == 1){
        if(p[[n]] == "NULL"){
          defaults[[n]] <- NULL
        }
      }
    }else{
      warning("set_parameters: unknown parameter", n, " in parameters file/object\n")
    }
  }
  a <- a[-grep("experiment", names(a))]
  a <- a[-grep("parameters", names(a))]
  for(n in names(a)){
    defaults[[n]] <- a[[n]]
  }
  defaults$dissimilarities <- as.character(defaults$dissimilarities)
  print(defaults)

  return(defaults)
}


#' check_end_generate_data
#'
#' @param experiment An mg-class object
#' @param analyze A vector of strings specifying which features should be analyzed
#' @param analyzes A vector of strings specifying which analyzes should be performed
#' @param dissimilarities A vector of strings specifying at which dissimilarities OTUs should be constructed
#'
#' @return An mg-class objects with OTUs/PICRUSt2 results
#' @export check_and_generate_data
#'
#'

check_and_generate_data <- function(
    experiment = NULL,
    analyze = NULL,
    analyzes = NULL,
    dissimilarities = NULL,
    processors=1){


  # This function checks if OTUs are to be analyzed and if so, whether they need to be constructed for given dissimilarities,
  # and if PICRUSt2 results are needed and if so, whether they are ready.
  result <- experiment
  basename <- experiment$metadata$basename
  ### construct OTUs if needed, only for missing dissimilarities (therefore one by one, and not in one go - to avoid overwriting those present in the object)
  if(sum(grepl('OTUs', analyze))){
    message("analyzeAll: checking if OTUs need to be constructed")
    todo <- c()
    for(diss in dissimilarities){
      if(is.null(experiment$OTUs$dists[[diss]])){
        message("analyzeAll: no OTUs for ", diss, ", they will be constructed")
        todo <- c(todo, diss)
      }
    }
    if(length(todo > 0)){
      # first check if we have working mothur
      mothur <- check_mothur(experiment=experiment)

      for(diss in todo){
        message("analyzeAll: generating OTUs for ", diss, " dissimilarities")
        if(!file.exists(experiment$ASVs$fasta.file)){
          write.table(experiment$ASVs$fasta, experiment$ASVs$fasta.file, quote=F, sep="", col.names=F)
          message("analyzeAll: fasta sequences written to ", experiment$ASVs$fasta.file)
        }
        if(!file.exists(experiment$ASVs$count.file)){
          write.table(experiment$ASVs$count_table, experiment$ASVs$count.file, quote=F, sep="\t", row.names=F)
          message("analyzeAll: count table written to ", experiment$ASVs$count.file)
        }
        message("analyzeAll: ", experiment$metadata$tree.algorithm, " will be used to construct a tree")
        result$OTUs$dists[[diss]] <- mgmisc1::constructOTUs(
          experiment = experiment,
          diss = diss,
          processors = processors
        )[[diss]]
        message("analyzeAll: OTUs generated")
        setwd(outdir)
        if(save){
          name <- paste0(experiment$metadata$basename, "_mg.rds")
          message("analyzeAll: saving experiment with OTUs as ", name)
          saveRDS(result, name)
        }
      }
    }else{
      message("analyzeAll: OTUs are present for all dissimilarities")
    }
  }

  # if KOs or pathways are to be analyzed they need to be present, check and generate if needed

  if(experiment$metadata$marker == "16S"){
    if(sum(grepl("KOs|pathways", analyze))){
      message("analyzeAll: checking if PICRUSt2 data are present")
      picrust.outdir <- paste0(experiment$metadata$basename, "_picrust2_outdir")
      picrust.filename.base <- paste0(experiment$metadata$basename, "_picrust2")

      if(length(mgmisc1::extractFeatureSet(experiment=experiment, feature="KOs", rarefied=F)) == 0 ||  length(mgmisc1::extractFeatureSet(experiment=experiment, feature="pathways", rarefied=F)) == 0){
        message("analyzeAll: no PICRUSt2 data, generating")
        # first check if we know how to run PICRUSt2
        picrust2 <- check_picrust(experiment=experiment)

        fs <- mgmisc1::extractFeatureSet(experiment=experiment, feature="ASVs", rarefied=F)
        if(processors <= 8){
          picrust.processors <- processors
        }else{
          picrust.processors <- 8
        }
        if(dir.exists(picrust.outdir)){
          message("analyzeAll: PICRUSt output directory found, trying to read in computation results")
          p <- mgmisc1::readPICRUSt2Dir(picrust.outdir)
          p$sampledata <- fs$sampledata
          if(!p$success){
            message("analyzeAll: ouput directory for PICRUSt2 does not contain full computation results, deleting")
            unlink(picrust.outdir, recursive=T)
            message("analyzeAll: computing PICRUSt2 predictions for non-rarefied data")
            p <- mgmisc1::preparePICRUSt2(featuretab=fs$featuretab,
                                          sampledata=fs$sampledata,
                                          fasta.file=fasta.file,
                                          filename.base=picrust.filename.base,
                                          conda=conda, condaenv=condaenv,
                                          picrust=picrust,
                                          picrust.outdir=picrust.outdir,
                                          plot.nsti=F,
                                          variable=NULL,
                                          processors=picrust.processors)
          }
          message("analyzeAll: PICRUSt2 predictions for non-rarefied data read in")
        }else{
          message("analyzeAll: computing PICRUSt2 predictions for non-rarefied data")
          p <- mgmisc1::preparePICRUSt2(featuretab=fs$featuretab,
                                        sampledata=fs$sampledata,
                                        fasta.file=fs$fasta.file,
                                        filename.base=picrust.filename.base,
                                        conda=picrust2$conda,
                                        condaenv=picrust2$condaenv,
                                        picrust=picrust2$picrust,
                                        picrust.outdir=picrust.outdir,
                                        plot.nsti=F,
                                        variable=NULL,
                                        processors=picrust.processors)
        }
        message("analyzeAll: number of rows in sampledata: ", nrow(p$sampledata))
        fname <- "picrust_nonrarefied_results"
        assign(fname, p, envir = parent.frame())
        result$KOs$featuretab <- p$KOs$featuretab
        result$KOs$featureannot <- p$KOs$featureannot
        result$KOs$sampledata <- p$sampledata
        result$KOs$nsti <- p$nsti
        result$pathways$featuretab <- p$pathways$featuretab
        result$pathways$featureannot <- p$pathways$featureannot
        result$pathways$sampledata <- p$sampledata
      }
    }
  }else{
    message("check_and_generate_data: PICRUSt2 can't be run for non-16S data")
  }

  # if alpha- or beta-diversity is to be analyzed, rarefied data are needed - check if present and generate if needed
  if( sum(grepl("alpha-diversity|beta-diversity", analyzes))){
    if(length(mgmisc1::extractFeatureSet(experiment=experiment, feature="ASVs", rarefied=T)) == 0){
      message("analyzeAll: no rarefied data for ASVs, generating")
      fasta.file <- sub("\\.fasta", ".rarefied.fasta", experiment$ASVs$fasta.file)
      tree.file <- sub("\\.tre", ".rarefied.tre", experiment$ASVs$tree.file)
      count.file <- sub("\\.count_table", ".rarefied.count_table", experiment$ASVs$count.file)
      result$ASVs$rarefied.data <- mgmisc1::rarefyData(featuretab=experiment$ASVs$featuretab,
                                                           featureannot=experiment$ASVs$featureannot,
                                                           sampledata=experiment$ASVs$sampledata,
                                                           tree=experiment$ASVs$tree,
                                                           tree.file=tree.file,
                                                           fasta=experiment$ASVs$fasta,
                                                           fasta.file=fasta.file,
                                                           count_table=experiment$ASVs$count_table,
                                                           count.file=count.file,
                                                           rarefy=1000, # no value set in experiment, therefore default
                                                           depth="auto")
    }
    for( feature in analyze ){
      message("check_and_generate_data: analyzing ", feature)
      if(feature == "OTUs"){
        for( diss in names(experiment$OTUs$dists) ){
          if(length(mgmisc1::extractFeatureSet(experiment=experiment, feature=feature, rarefied=T, d=diss)) == 0){
            message("analyzeAll: no rarefied data for OTUs at ", diss, " dissimilarity, generating")
            fs <- mgmisc1::extractFeatureSet(experiment=experiment, feature="OTUs", rarefied=F, d=diss)
            fasta.file <- sub("\\.fasta", ".rarefied.fasta", fs$fasta.file)
            tree.file <- sub("\\.tre", ".rarefied.tre", fs$tree.file)
            count.file <- sub("\\.count_table", ".rarefied.count_table", fs$count.file)
            result$OTUs$dists[[diss]]$rarefied.data <- mgmisc1::rarefyData(featuretab=fs$featuretab,
                                                                      featureannot=fs$featureannot,
                                                                      sampledata=fs$sampledata,
                                                                      tree=fs$tree,
                                                                      tree.file=tree.file,
                                                                      fasta=fs$fasta,
                                                                      fasta.file=fasta.file,
                                                                      count_table=fs$count_table,
                                                                      count.file=count.file,
                                                                      rarefy=experiment$ASVs$rarefied.data$rarefy,
                                                                      depth=experiment$ASVs$rarefied.data$depth)
          }
        }
      }else{
        if(length(mgmisc1::extractFeatureSet(experiment=experiment, feature="KOs", rarefied=T)) == 0 || length(mgmisc1::extractFeatureSet(experiment=experiment, feature="pathways", rarefied=T)) == 0 ){
          if(experiment$metadata$marker == "16S"){
          message("check_and_generate_data: no rarefied PICRUSt2 data, generating")
          picrust2 <- check_picrust(experiment=experiment)
          fs <- mgmisc1::extractFeatureSet(experiment=experiment, feature="ASVs", rarefied=T)
          picrust.outdir <- paste0(basename, "_rarefied_picrust2_outdir")
          picrust.filename.base <- paste0(basename, "_rarefied_picrust2")
          if(processors <= 8){
            picrust.processors <- processors
          }else{
            picrust.processors <- 8
          }
          if(dir.exists(picrust.outdir)){
            message("check_and_generate_data: PICRUSt output directory found, trying to read in computation results")
            p <- mgmisc1::readPICRUSt2Dir(picrust.outdir)
            p$sampledata <- fs$sampledata
            if(!p$success){
              message("check_and_generate_data: ouput directory for PICRUSt2 does not contain full computation results, deleting")
              unlink(picrust.outdir, recursive=T)
              message("check_and_generate_data: computing PICRUSt2 predictions for rarefied data")
              p <- mgmisc1::preparePICRUSt2(featuretab=fs$featuretab,
                                            sampledata=fs$sampledata,
                                            fasta.file=fs$fasta.file,
                                            filename.base=picrust.filename.base,
                                            conda=picrust2$conda,
                                            condaenv=picrust2$condaenv,
                                            picrust=picrust2$picrust,
                                            picrust.outdir=picrust.outdir,
                                            plot.nsti=F,
                                            variable=NULL,
                                            processors=picrust.processors)
            }
            message("check_and_generate_data: PICRUSt2 predictions for rarefied data read in")
          }else{
            message("check_and_generate_data: computing PICRUSt2 predictions for rarefied data")
            p <- mgmisc1::preparePICRUSt2(featuretab=fs$featuretab,
                                          sampledata=fs$sampledata,
                                          fasta.file=fs$fasta.file,
                                          filename.base=picrust.filename.base,
                                          conda=picrust2$conda,
                                          condaenv=picrust2$condaenv,
                                          picrust=picrust2$picrust,
                                          picrust.outdir=picrust.outdir,
                                          plot.nsti=F,
                                          variable=NULL,
                                          processors=picrust.processors)
          }
          fname <- "picrust_rarefied_results"
          assign(fname, p, envir = parent.frame())
          message("check_and_generate_data: number of rows in sampledata: ", nrow(p$sampledata))
          message("check_and_generate_data: placing PICRUSt2 results in the experiment object")
          result$KOs$rarefied.data$featuretab <- p$KOs$featuretab
          message("featuretab")
          result$KOs$rarefied.data$featureannot <- p$KOs$featureannot
          message("featureannot")
          result$KOs$rarefied.data$sampledata <- p$sampledata
          message("sampledata")
          result$KOs$rarefied.data$nsti <- p$nsti
          message("nsti")
          result$pathways$rarefied.data$featuretab <- p$pathways$featuretab
          message("featuretab")
          result$pathways$rarefied.data$featureannot <- p$pathways$featureannot
          message("featureannot")
          result$pathways$rarefied.data$sampledata <- p$sampledata
          message("sampledata")
        }
        message("1")
        }
      }
      message("2")
    }
    message("3")
  }

  message("check_and_generate_data: all results placed in the experiment object, returning the object")
  return(result)
}



check_picrust <- function(experiment=NULL){
  # the function returns a list with conda.path, condaenv, and picrust elements being respective paths to conda and
  # picrust2_pipeline.py as well as the name of conda environment in which PICRUSt2 is installed, plus conda and PICRUSt2 versions
  res <- list()

  res$conda.path <- ""
  res$conda.version <- ""
  res$condaenv <- ""
  res$picrust <- ""
  res$picrust.version <- ""

  # is there picrust2_pipeline.py in experiment's metadata?
  message("analyzeAll: checking for PICRUSt2 installation...")
  if(!((is.null(experiment$metadata$conda) || experiment$metadata$conda == "") || (is.null(experiment$metadata$condaenv) || experiment$metadata$condaenv == "") || (is.null(experiment$metadata$picrust) || experiment$metadata$picrust == ""))){
    # there are paths in experiment's metadata, check if they work
    # is conda working?
    conda.version <- try(system2(experiment$metadata$conda, "--version", stdout=T))
    if(!inherits(conda.version, "try-error")){
      # conda works
      message("analyzeAll: working conda found in experiment's metadata")
      res$conda.path <- experiment$metadata$conda
      res$conda.version <- conda.version
    }else{
      # conda doesn't work, find working binary
      message("analyzeAll: no working conda in experiment's metadata, searching elsewhere")
      conda <- find_conda()
      # if we get that far, we have working conda
      message("analyzeAll: working conda v.", conda$conda.version, " found in ", conda$conda.path)
      res$conda.path <- conda$conda.path
      res$conda.version <- conda$conda.version
    }

    # we have working conda, check if PICRUSt2 is installed in the environment given in metadata
    message("analyzeAll: checking if PICRUSt2 is installed in conda environment give in experiment's metadata")
    if(reticulate::condaenv_exists(experiment$metadata$condaenv)){
      # there is an environment given in metadata, check if it contains PICRUSt2 installation
      reticulate::use_condaenv(conda=res$conda.path, condaenv=experiment$metadata$condaenv, required=T)
      reticulate::py_config() # python initialization, setting up env vars so that contents of conda environment's bin directory could be accessed
      picrust.version <- try(system2(experiment$metadata$picrust, "--version", stdout=T))
      if(!inherits(picrust.version, "try-error")){
        # PICRUSt2 works, all is fine
        message("analyzeAll: working PICRUSt2 found in \"", experiment$metadata$condaenv, "\"")
        res$condaenv <- experiment$metadata$condaenv
        res$picrust <- experiment$metadata$picrust
        res$picrust.version <- picrust.version
        return(res)
      }else{
        # no working PICRUSt in a given environment, try finding another env with PICRUSt2
        message("analyzeAll: there is no working PICRUSt2 installation in conda environment \"", experiment$metadata$condaenv, "\", need to search in other envs")
        picrust.env <- find_picrust(conda=res$conda.path, picrust=experiment$metadata$picrust)
        if(length(picrust.env) != 0){
          res$condaenv <- picrust.env$condaenv
          res$picrust <- picrust.env$picrust
          res$picrust.version <- picrust.env$picrust.version
        }
      }
    }else{
      # there is no such an environment, check if there is any env with PICRUSt2 installed (assume standard name of PICRUSt2 script - 'picrust2_pipeline.py')
      message("analyzeAll: there is no conda environment named \"", experiment$metadata$condaenv, "\", need to search for PICRUSt2 in other envs")
      picrust.env <- find_picrust(conda=res$conda.path, picrust="")
      if(length(picrust.env) != 0){
        res$condaenv <- picrust.env$condaenv
        res$picrust <- picrust.env$picrust
        res$picrust.version <- picrust.env$picrust.version
      }
    }
  }else{
    # at least one path is missing, try to find conda and environment in which PICRUSt2 is installed
    # find working conda
    message("analyzeAll: some paths in experiment's metadata are missing, looking for conda")
    conda <- find_conda()
    message("analyzeAll: working conda v.", conda$conda.version, " found in ", conda$conda.path)
    res$conda.path <- conda$conda.path
    res$conda.version <- conda$conda.version

    # find an environment in which PICRUSt2 is installed
    message("analyzeAll: looking for PICRUSt2 installation")
    picrust.env <- find_picrust(conda=res$conda.path, picrust="")
    if(length(picrust.env) != 0){
      res$condaenv <- picrust.env$condaenv
      res$picrust <- picrust.env$picrust
      res$picrust.version <- picrust.env$picrust.version
    }
  }

  return(res)
}

find_conda <- function(){
  res <- list()
  conda.path <- try(reticulate::conda_binary())
  if(inherits(conda.path, "try-error")){
    # no conda can be found, ask for a path, maybe it lives in a non-standard location
    message("analyzeAll: no working conda binary could be found")
    while(TRUE){
      cat("Type in a valid path to conda binary: ")
      conda.path <- readLines(n=1)
      conda.version <- try(system2(conda.path, "--version", stdout=T))
      if(inherits(conda.version, "try-error")){
        message(conda.path, " is not a path to working conda binary")
      }else{
        # we have a path to working conda binary
        message(conda.path, "is a working conda binary")
        break
      }
    }
    res$conda.path <- conda.path
    res$conda.version <- conda.version
  }else{
    conda.version <- try(system2(conda.path, "--version", stdout=T))
    message("analyzeAll: working conda v. ", conda.version, " found in ", conda.path)
    res$conda.path <- conda.path
    res$conda.version <- conda.version
  }
  return(res)
}

find_picrust <- function(conda, picrust){
  # need to add code to compare versions and choose the latest one
  res <- list()
  if(picrust == ""){
    picrust <- "picrust2_pipeline.py"
  }

  #is picrust2_pipeline.py in PATH? A suitable environment might have been already activated
  program_path <- system2("which", picrust, stdout=T)
  if(length(program_path) != 0){
    # picrust2_pipeline.py is in PATH, check if it is in some conda env
    env.name <- NULL
    if(grepl("/envs/", program_path)){
      # picrust2 script is installed in an active conda env, we need to find the env's name
      env.name <- sub("/bin.*", "", sub(".*/envs/", "", program_path))
    }
    res$condaenv <- env.name # "" if picrust2_pipline.py is not in any conda env
    res$picrust <- picrust
    res$picrust.version <- system2(program_path, "--version", stdout=T)
    return(res)
  }

  # get a list of all environments
  all_envs <- reticulate::conda_list(conda=conda)

  # for each env check if PICRUSt2 is installed
  for(env.name in all_envs$name){
    message("find_picrust: checking ", env.name, " environment")
    env.bin.path <- sub("/python$", "", all_envs[ all_envs$name == env.name, ]$python)
    print(env.bin.path)
    program_path <- try(system2("which", picrust, stdout = TRUE, stderr = FALSE, env = paste("PATH=", env.bin.path, sep="")), silent = TRUE)
    if(!inherits(program_path, "try-error") && length(program_path) != 0){
      message("find_picrust: found PICRUSt2 in ", env.name)
      res$condaenv <- env.name
      res$picrust <- picrust
      res$picrust.version <- system2(program_path, "--version", stdout=T)
      break
    }
  }

  return(res)
}
