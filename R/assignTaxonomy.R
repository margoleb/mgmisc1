function (seqs,
          refFasta,
          minBoot = 50,
          tryRC = FALSE,
          outputBootstraps = FALSE,
          taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family",
                        "Genus", "Species"),
          multithread = FALSE,
          verbose = FALSE,
          depth = 3)
{
  MIN_REF_LEN <- 20
  MIN_TAX_LEN <- 50
  seqs <- getSequences(seqs)
  if (min(nchar(seqs)) < MIN_TAX_LEN) {
    warning("Some sequences were shorter than ", MIN_TAX_LEN,
            " nts and will not receive a taxonomic classification.")
  }
  refsr <- readFasta(refFasta)
  lens <- width(sread(refsr))
  if (any(lens < MIN_REF_LEN)) {
    refsr <- refsr[lens >= MIN_REF_LEN]
    warning(paste0("Some reference sequences were too short (<",
                   MIN_REF_LEN, "nts) and were excluded."))
  }
  refs <- as.character(sread(refsr))
  tax <- as.character(id(refsr))
  tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x))
  UNITE <- FALSE
  if (all(grepl("FU\\|re[pf]s", tax[1:10]))) {
    UNITE <- TRUE
    cat("UNITE fungal taxonomic reference detected.\n")
    tax <- sapply(strsplit(tax, "\\|"), `[`, 5)
    tax <- gsub("[pcofg]__unidentified;", "_DADA2_UNSPECIFIED;",
                tax)
    tax <- gsub(";s__(\\w+)_", ";s__", tax)
    tax <- gsub(";s__sp$", ";_DADA2_UNSPECIFIED", tax)
  }
  if (!grepl(";", tax[[1]])) {
    if (length(unlist(strsplit(tax[[1]], "\\s"))) == 3) {
      stop("Incorrect reference file format for assignTaxonomy (this looks like a file formatted for assignSpecies).")
    }
    else {
      stop("Incorrect reference file format for assignTaxonomy.")
    }
  }
  #tax.depth <- sapply(strsplit(tax, ";"), length)
  td <- depth
  for (i in seq(length(tax))) {
    if (tax.depth[[i]] < td) {
      for (j in seq(td - tax.depth[[i]])) {
        tax[[i]] <- paste0(tax[[i]], "_DADA2_UNSPECIFIED;")
      }
    }
  }
  genus.unq <- unique(tax)
  ref.to.genus <- match(tax, genus.unq)
  tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol = td,
                    byrow = TRUE)
  tax.df <- as.data.frame(tax.mat)
  for (i in seq(ncol(tax.df))) {
    tax.df[, i] <- factor(tax.df[, i])
    tax.df[, i] <- as.integer(tax.df[, i])
  }
  tax.mat.int <- as.matrix(tax.df)
  if (is.logical(multithread)) {
    if (multithread == TRUE) {
      RcppParallel::setThreadOptions(numThreads = "auto")
    }
    else {
      RcppParallel::setThreadOptions(numThreads = 1)
    }
  }
  else if (is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
  }
  else {
    warning("Invalid multithread parameter. Running as a single thread.")
    RcppParallel::setThreadOptions(numThreads = 1)
  }
  assignment <- C_assign_taxonomy2(seqs, rc(seqs), refs, ref.to.genus,
                                   tax.mat.int, tryRC, verbose)
  bestHit <- genus.unq[assignment$tax]
  boots <- assignment$boot
  taxes <- strsplit(bestHit, ";")
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,
  ] >= minBoot])
  tax.out <- matrix(NA_character_, nrow = length(seqs), ncol = td)
  for (i in seq(length(seqs))) {
    if (length(taxes[[i]]) > 0) {
      tax.out[i, 1:length(taxes[[i]])] <- taxes[[i]]
    }
  }
  rownames(tax.out) <- seqs
  colnames(tax.out) <- taxLevels[1:ncol(tax.out)]
  tax.out[tax.out == "_DADA2_UNSPECIFIED"] <- NA_character_
  if (outputBootstraps) {
    boots.out <- matrix(boots, nrow = length(seqs), ncol = td)
    rownames(boots.out) <- seqs
    colnames(boots.out) <- taxLevels[1:ncol(boots.out)]
    list(tax = tax.out, boot = boots.out)
  }
  else {
    tax.out
  }
}
