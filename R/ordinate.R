#' @export ordinate
#' @import vegan

ordinate <- function(
  method = NULL,
  d = NULL,
  distname = NULL,
  featuretab = NULL,
  featureannot = NULL,
  sampledata = NULL,
  shape.variable = NULL,
  sample.size.variable = NULL,
  biplot = F,
  feature.color.variable = NULL,
  size.means.abundance = F,
  abundance.cutoff = 50,
  triplot = F,
  chemistry.variables = NULL,
  sqrt.dist = F,
  statistics = F,
  varpart.variable = NULL,
  rand.seed = 667,
  p.value = 0.05,
  permu = 999, # number of permutations for tests
  steps = 200, # max number of steps in ordistep
  ordistep.direction = 'both' # forward, reverse or both types of selection in ordistep?
){
  methods <- c('PCA', 'MDS', 'NMDS', 'RDA', 'dbRDA', 'CCA')
  if(! method %in% methods){
    stop("ordinate: unknown type of ordination\n method must be one of: PCA, MDS, NMDS, RDA, dbRDA or CCA")
  }

  dfsamples <- NULL
  dffeatures <- NULL
  dfvectors <- NULL
  title=""

  if(method %in% c('RDA', 'dbRDA', 'CCA')){
    if(triplot){
      chem <- sampledata[ , colnames(sampledata) %in% chemistry.variables ]
      chem$shape <- sampledata[[shape.variable]]
      # 'shotgun' model selection with ordistep
      if( method == 'dbRDA'){
        set.seed(rand.seed)
        mod0 <- vegan::dbrda(d ~ 1, data=chem, sqrt.dist=sqrt.dist)
        mod1 <- vegan::dbrda(d ~ ., data=chem, sqrt.dist=sqrt.dist)
      }else if(method == 'RDA'){
        set.seed(rand.seed)
        mod0 <- vegan::rda(featuretab ~ 1, data=chem)
        mod1 <- vegan::rda(featuretab ~ ., data=chem)
      }else{
        set.seed(rand.seed)
        mod0 <- vegan::cca(featuretab ~ 1, data=chem)
        mod1 <- vegan::cca(featuretab ~ ., data=chem)
      }
      set.seed(rand.seed)
      ord <- ordistep(mod0, scope=formula(mod1), direction=ordistep.direction, permu=permu, steps=steps, Pin=p.value, Pout=2*p.value)
      dfvectors <- as.data.frame(scores(ord, display='bp'))
      colnames(dfvectors) <- c('x', 'y')
      dfvectors <- dfvectors[ rownames(dfvectors) %in% chemistry.variables, ]
    }else{ # no triplot to be generated
      if(biplot){
        if(method == 'dbRDA'){
          set.seed(rand.seed)
          f <- as.formula(paste0("d ~ ", shape.variable))
          ord <- try(vegan::dbrda( f, data=sampledata, sqrt.dist=sqrt.dist ), silent=T)
          if(inherits(ord, "try-error")) {
            set.seed(rand.seed)
            ord <- vegan::capscale(f, data=sampledata)
            vegan::sppscores(ord) <- sqrt(vegan::decostand(featuretab, "total"))
            outlist$dbrda <- ord
            dffeatures <- as.data.frame(vegan::scores(ord, display='species'))
          }
        }else if(method == 'RDA'){
          set.seed(rand.seed)
          ord <- vegan::rda(featuretab ~ shape.variable, data=sampledata)
        }else{
          set.seed(rand.seed)
         ord <- vegan::cca(featuretab ~ shape.variable, data=sampledata)
        }
      }
    }
  }else{
    if(method == 'PCA'){ # principal component analysis - RDA with no constraints
      set.seed(rand.seed)
      ord <- vegan::rda(featuretab ~ 1, data=sampledata)
      if(biplot){

      }
    }else if(method == 'MDS'){ # a.k.a. 'PCoA' principal coordinates analysis - dbRDA with no constraints
      set.seed(rand.seed)
      ord <- vegan::dbrda(d ~ 1, data=sampledata, sqrt.dist=sqrt.dist)
      if(biplot){

      }
    }else if(method == 'NMDS'){ # NMDS
      set.seed(rand.seed)
      ord <- vegan::metaMDS(d, k=2, try=try, trymax=trymax)
      if(biplot){
        ord.points <- postMDS(vegan::scores(ord, display='sites', choices=c(1,2)), d)
        dffeatures <- as.data.frame(wascores(ord.points, featuretab))
      }
    }else{
      warning("ordinate: unknown ordination method ", method, ", skipping")
    }
  }

  outlist$ord <- ord
  outlist$shape.variable <- shape.variable
  outlist$sample.color.variable <- sample.color.variable
  outlist$varpart.variable <- varpart.variable
  outlist$rand.seed <- rand.seed
  outlist$varpart <- NULL

  dfsamples <- as.data.frame(vegan::scores(ord, display='sites'))
  colnames(dfsamples) <- c('x', 'y')
  dfsamples$shape <- sampledata[[shape.variable]]
  dfsamples$color <- sampledata[[sample.color.variable]]
  if(!is.null(sample.size.variable)){
    dfsamples$size <- sampledata[[sample.size.variable]]
    if(!is.numeric(dfsamples$size)){
      dfsamples$size <- as.numeric(levels(as.factor(dfsamples$size)))
    }
  }

  if(!is.null(varpart.variable)){
    message("ordinate: Performing variance partitioning")
    if(triplot){
      chem.vars <- paste(chemistry.variables, collapse=" + ")
      f1 <- as.formula(paste0(" ~ ", shape.variable))
      f2 <- as.formula(paste0(" ~ ", chem.vars))
      ordvarpart <- vegan::varpart(d, f1, f2, data=sampledata)
      explvar.shape <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[1,3])), '%');
      title <- paste0("\nVariance explained by ", shape.variable, ": ", explvar.shape );
      explvar.chem <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[2,3])), '%')
      outlist$varpart <- ordvarpart
      title <- paste0(title, "\nVariance explained by physicochemical variables: ", explvar.chem );
    }else{
      f1 <- as.formula(paste0(" ~ ", shape.variable))
      f2 <- as.formula(paste0(" ~ ", varpart.variable))
      ordvarpart <- vegan::varpart(d, f1, f2, data=sampledata)
      explvar <- paste0(sprintf("%.2f", 100 * (ordvarpart$part$indfract[1,3])), '%');
      title <- paste0(title, "\nVariance explained by ", shape.variable, ": ", explvar );
      outlist$varpart <- ordvarpart
    }
  }

  if(biplot){
    colnames(dffeatures) <- c('x', 'y')
    dffeatures$color <- featureannot[[feature.color.variable]]
    dffeatures$Abundance <- colSums(featuretab)
    dffeatures <- head(dffeatures[ order(dffeatures$Abundance, decreasing=T), ], abundance.cutoff)
  }

  outlist$betadisper <- NULL
  outlist$betadisper.anova <- NULL
  outlist$permutest <- NULL

  if( as.logical(statistics) ) {
    message("ordinate: Testing model and performing betadisper analysis")
    ordbetadisper <- vegan::betadisper(as.dist(d), sampledata[[shape.variable]], bias.adjust=T) # in case of GUniFrac-produced matrices they need to be transformed to proper dist-class objects
    outlist$betadisper <- ordbetadisper
    ordbetadisperanova <- anova(ordbetadisper, permu=permu)
    outlist$betadisper.anova <- ordbetadisperanova
    betadisperf <- sprintf( "%.2f", ordbetadisperanova$"F value"[1] );
    betadisperp <- sprintf( "%.2e", ordbetadisperanova$"Pr(>F)"[1] );
    if(method %in% c('CCA', 'dbRDA', 'RDA')){
      orddbrdaanova <- vegan::anova(ord, permu=permu)
      outlist$permutest <- orddbrdaanova
      dbrdaanovaf <- sprintf( "%.2f",orddbrdaanova$F[1] );
      dbrdaanovap <- sprintf( "%.3f", orddbrdaanova$"Pr(>F)"[1] );
      title <- paste0("ANOVA F = ", dbrdaanovaf, ", p = ", dbrdaanovap, "\nbetadisper F = ", betadisperf, ", p = ", betadisperp, title)
    }else{
      f <- as.formula(paste0("d ~ ", shape.variable))
      permutest <- vegan::adonis2(f, data=sampledata, permu=permu)
      outlist$permutest <- permutest
      message("plotNMDS: PERMANOVA done")
      permanovaf <- sprintf( "%.2f", permutest$F[1] );
      permanovap <- sprintf( "%.3f", permutest$"Pr(>F)"[1] );
      title <- paste0("ANOVA F = ", permanovaf, ", p = ", permanovap, "\nbetadisper F = ", betadisperf, ", p = ", betadisperp, title)
    }
  }


  outlist <- list(ord=ord, permutest=permutest, betadisper=betadisper, betadisper.anova=betadisper.anova, dfsamples=dfsamples, dffeatures=dffeatures, dfvectors=dfvectors, varpart=ordvarpart, title=title)
  return(outlist)
}
