predictOneGene <- function(gene, bamfiles, fitpar, genome=Hsapiens,
                           models, readlength, minsize, maxsize) {
  stopifnot(is(gene, "GRanges"))
  stopifnot(all(sapply(models, function(x) names(x) %in% c("formula","offset"))))
  stopifnot(!is.null(fitpar))
  stopifnot(all(names(bamfiles) %in% names(fitpar)))
  if (is.null(names(bamfiles))) {
    names(bamfiles) <- seq_along(bamfiles)
  }
  fragtypes <- buildFragtypesFromExons(gene, genome, readlength=readlength,
                                       minsize=minsize, maxsize=maxsize)
  res <- list()
  for (bamname in names(bamfiles)) {
    # add counts
    bamfile <- bamfiles[bamname]
    generange <- range(gene)
    strand(generange) <- "*" # not necessary
    suppressWarnings({
      ga <- readGAlignAlpine(bamfile, generange)
    })
    if (length(ga) == 0) {
      res[[bamname]] <- as.list(rep(NA,length(models)))
      names(res[[bamname]]) <- names(models)
      next
    }
    ga <- keepSeqlevels(ga, as.character(seqnames(gene)[1]))
    fco <- findCompatibleOverlaps(ga, GRangesList(gene))
    # message("-- ",round(length(fco)/length(ga),2)," compatible overlaps")
    reads <- gaToReadsOnTx(ga, GRangesList(gene), fco)

    # save fragment coverage for later
    l <- sum(width(gene))
    frag.cov <- coverage(reads[[1]][start(reads[[1]]) != 1 & end(reads[[1]]) != l])
    
    fragtypes.temp <- matchReadsToFraglist(reads, list(fragtypes))[[1]]
    ## -- fragment bias --
    fraglen.density <- fitpar[[bamname]][["fraglen.density"]]
    stopifnot(!is.null(fraglen.density))
    fragtypes.temp$logdfraglen <- log(matchToDensity(fragtypes.temp$fraglen, fraglen.density))
    ## -- random hexamer priming bias with VLMM --
    vlmm.fivep <- fitpar[[bamname]][["vlmm.fivep"]]
    vlmm.threep <- fitpar[[bamname]][["vlmm.threep"]]
    stopifnot(!is.null(vlmm.fivep))
    stopifnot(!is.null(vlmm.threep))
    fragtypes.temp <- addVLMMBias(fragtypes.temp, vlmm.fivep, vlmm.threep)
    
    # -- fit models --
    res[[bamname]] <- list()

    # remove first and last bp for predicting coverage along transcript
    not.first.or.last.bp <- !(fragtypes.temp$start == 1 | fragtypes.temp$end == l)
    fragtypes.temp <- fragtypes.temp[not.first.or.last.bp,]

    # summarize counts and predicted counts on the fragment midpoint
    ir <- IRanges(fragtypes.temp$start, fragtypes.temp$end)
    fac <- factor(fragtypes.temp$mid, seq_len(l))
    cts <- sapply(split(fragtypes.temp$count, fac), sum)
    gc <- sapply(split(fragtypes.temp$gc, fac), mean)
    res[[bamname]]$l <- l
    res[[bamname]]$frag.cov <- frag.cov
    res[[bamname]]$cts <- cts
    res[[bamname]]$gc <- gc
    res[[bamname]]$pred.cts <- list()
    res[[bamname]]$pred.cov <- list()
    for (modeltype in names(models)) {
      # message("predicting model type: ",modeltype)
      log.lambda <- getLogLambda(fragtypes.temp, models, modeltype, fitpar, bamname)
      pred0 <- exp(log.lambda)
      pred <- pred0/mean(pred0)*mean(fragtypes.temp$count)
      pred.cts <- sapply(split(pred, fac), sum)
      res[[bamname]][["pred.cts"]][[modeltype]] <- pred.cts
      res[[bamname]][["pred.cov"]][[modeltype]] <- coverage(ir, weight=pred)
    }
  }
  res
}
