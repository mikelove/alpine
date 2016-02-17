splitGenesAcrossChroms <- function(ebg, txdf) {
  split.chroms <- sapply(split(txdf$TXCHROM, txdf$GENEID), function(x) !all(x == x[1]))
  message("found ",sum(split.chroms),
          " genes split over chroms out of ",length(split.chroms))
  split.chroms <- names(split.chroms)[split.chroms]
  new.genes <- GRangesList()
  for (gid in split.chroms) {
    chroms <- unique(txdf$TXCHROM[txdf$GENEID == gid])
    exs <- ebg[[gid]]
    for (i in seq_along(chroms)) {
      # cl = chromosome split
      new.name <- paste0(gid,"_cs",i)
      txdf$GENEID[txdf$GENEID == gid & txdf$TXCHROM == chroms[i]] <- new.name
      new.genes[[new.name]] <- exs[seqnames(exs) == chroms[i]]
    }
    ebg[[gid]] <- NULL
  }
  ebg <- c(ebg, new.genes)
  list(ebg=ebg, txdf=txdf)
}
splitLongGenes <- function(ebg, ebt, txdf, long=1e6) {
  strand(ebg) <- "*"
  r <- unlist(range(ebg))
  w <- width(r)
  stopifnot(length(w) == length(ebg))
  long.genes <- names(ebg)[w > long]
  message("found ",length(long.genes),
          " long genes (1e",log10(long)," bp) out of ",length(ebg))
  new.genes <- GRangesList()
  for (gid in long.genes) {
    # ls = long split
    new.names <- paste0(gid,"_ls",seq_len(sum(txdf$GENEID == gid)))
    ebg[[gid]] <- NULL
    txdf$GENEID[txdf$GENEID == gid] <- new.names
    for (new.gene in new.names) {
      gr <- ebt[[txdf$TXNAME[txdf$GENEID == new.gene]]]
      mcols(gr)$exon_rank <- NULL
      new.genes[[new.gene]] <- gr
    }
  }
  ebg <- c(ebg, new.genes)
  list(ebg=ebg, txdf=txdf)
}
mergeGenes <- function(ebg, txdf) {
  fo <- findOverlaps(ebg, ignore.strand=TRUE)
  fo <- fo[queryHits(fo) < subjectHits(fo)]
  mat <- as.matrix(fo)
  graph <- ftM2graphNEL(mat, edgemode="undirected")
  components <- connectedComp(graph)
  message("found ",length(components), " clusters from ",length(ebg)," genes")
  components <- lapply(components, function(x) names(ebg)[as.numeric(x)])
  for (cluster in components) {
    txdf$GENEID[txdf$GENEID %in% cluster] <- paste0(cluster[1],"_mrg")
  }
  txdf
}
extractRes <- function(res, model, what, nsamp) {
  do.call(rbind, lapply(res, function(x) {
    if (is.null(x)) {
      if (what == "count") return(rep(0, nsamp))
      return(rep(NA, nsamp)) # the whole gene gets a single row of NA
    }
    if (what == "count") {
      return(sapply(x, `[[`, what))
    } else {
      res.list <- lapply(x, function(y) y[[model]][[what]])
      return(do.call(cbind, res.list))
    }
  }))
}
extract <- function(res, model, nsamp, lib.sizes=1e6) {
  fpkm <- extractRes(res, model, "theta", nsamp)
  lambda <- extractRes(res, model, "lambda", nsamp)
  count <- extractRes(res, model, "count", nsamp)
  lambdaBar <- colMeans(lambda, na.rm=TRUE)
  colSumsCount <- colSums(count)
  multFactor <-  lambdaBar * lib.sizes / colSumsCount
  sweep(fpkm, 2, multFactor, `*`)
}

###

nogenbank <- function(x) sub("(.*)\\(GenBank)","\\1",x)
getYmax <- function(res) {
  max(sapply(seq_along(res), function(i) max(res[[i]]$frag.cov)))
}
plotCovOne <- function(res, i=1, m="GC", xlab="", ylab="", log=FALSE, ...) {
  transform <- if (log) {
    function(x) log10(x + 1)
  } else {
    I
  }
  ymax <- transform(getYmax(res[i]))
  plot(transform(as.numeric(res[[i]]$frag.cov)), type="l", ylim=c(0,ymax),
       xlab=xlab, ylab=ylab, col=cond[i], lwd=3, ...)
  lines(transform(as.numeric(res[[i]][["pred.cov"]][[m]])), col=rgb(0,0,0,.7),lwd=3)
}
plotCov <- function(res, m="GC", cond, xlab="", ylab="", log=FALSE, lwd=3, ...) {
  for (i in seq_along(res)) {
    transform <- if (log) {
      function(x) log10(x + 1)
    } else {
      I
    }
    #ymax <- transform(getYmax(res))
    plot(transform(as.numeric(res[[i]]$frag.cov)), type="l", # ylim=c(0,ymax),
         xlab=xlab, ylab=ylab, col=cond[i], lwd=lwd, ...)
    lines(transform(as.numeric(res[[i]][["pred.cov"]][[m]])), col=rgb(0,0,0,.7),lwd=lwd)
  }
}
varExplainedPos <- function(res, m="GC", errorfun) {
  do.call(rbind,lapply(seq_along(res), function(i) {
    y <- as.numeric(res[[i]]$frag.cov)
    fit <- as.numeric(res[[i]][["pred.cov"]][[m]])
    y <- c(y, rep(0, res[[i]]$l - length(y)))
    fit <- c(fit, rep(0, res[[i]]$l - length(fit)))
    mu <- mean(y)
    null <- errorfun(y, mu)
    residual <- errorfun(y, fit)
    reduction <- if (null == 0) NA else (null - residual)/null
    data.frame(mu=mu, null=null, reduction=reduction)
  }))
}
getCountMatrix <- function(gene, bamfile, genome=Hsapiens) {
  fragtypes <- buildFragtypesFromExons(gene, genome)
  l <- sum(width(gene))
  generange <- range(gene)
  strand(generange) <- "*" # not necessary
  if (!as.character(seqnames(generange)) %in% seqlevels(BamFile(bamfile))) next
  flag <- alpineFlag()
  suppressWarnings({ ga <- readGAlignmentPairs(bamfile, param=ScanBamParam(which=generange,flag=flag)) })
  ga <- keepSeqlevels(ga, as.character(seqnames(gene)[1]))
  fco <- findCompatibleOverlaps(ga, GRangesList(gene))
  reads <- gaToReadsOnTx(ga, GRangesList(gene), fco)
  fraglist <- matchReadsToFraglist(reads, list(fragtypes))
  fragtypes <- fraglist[[1]]
  sp <- split(fragtypes$count, fragtypes$start)
  lens <- sapply(sp, length)
  maxl <- max(lens)
  short.idx <- unname(which(lens < maxl))
  for (i in short.idx) {
    sp[[i]] <- c(sp[[i]], rep(0, maxl - lens[i]))
  }
  mat <- do.call(cbind, sp)
  mat[nrow(mat):1,]
}
getFragStarts <- function(gene, bamfile, ends=FALSE, width=FALSE) {
  stopifnot(is(gene, "GRanges"))
  gene.grl <- GRangesList(gene)
  genestrand <- as.character(strand(gene[1]))
  generange <- range(gene)
  l <- sum(width(gene))
  strand(generange) <- "*" # not necessary
  if (!as.character(seqnames(generange)) %in% seqlevels(BamFile(bamfile))) return(NULL)
  flag <- alpineFlag()
  suppressWarnings({ ga <- readGAlignmentPairs(bamfile, param=ScanBamParam(which=generange,flag=flag)) })
  ga <- keepSeqlevels(ga, as.character(seqnames(gene)[1]))
  if (length(ga) == 0) return(NULL)
  reads <- gaToReadsOnTx(ga, gene.grl)[[1]]
  if (width) {
    return(width(reads))
  } else if (ends) {
    return(end(reads))
  } else {
    return(start(reads))
  }
}
plotFragStarts <- function(gene, bamfile, genome, window=200, ends=FALSE, ...) {
  x <- getFragStarts(gene, bamfile, ends=ends)
  if (is.null(x)) return(NULL)
  l <- sum(width(gene))
  reads <- as.numeric(table(factor(x, 1:l)))
  plot(reads, xlab="position", type="h", ...)
  if (!missing(genome)) {
    exon.dna <- unlist(getSeq(genome, gene))
    gc <- letterFrequencyInSlidingView(exon.dna, view.width=window, letters="GC", as.prob=TRUE)
    lines(seq_along(gc) + window/2, gc * max(reads), col="dodgerblue", lwd=2)
    axis(4, 0:4/4 * max(reads), 0:4/4, col="dodgerblue")
  }
}
roughSummarizeOverlaps <- function(features, bamfile) {
  # assume paired end
  rngs <- unlist(range(features))
  rngs$id <- paste(seqnames(rngs), start(rngs), end(rngs), sep="-")
  data <- countBam(bamfile, param=ScanBamParam(which=rngs))
  data$id <- paste(data$space, data$start, data$end, sep="-")
  data <- data[match(rngs$id, data$id),]
  round(data$records/2)
}
countBamfiles <- function(bamfiles, which) {
  # assume paired end
  round(sapply(seq_along(bamfiles), function(i) countBam(bamfiles[i], param=ScanBamParam(which=which))$records/2))
}
getReadlength <- function(bamfiles) {
  getRL1 <- function(file) {
    qwidth(readGAlignments(BamFile(file, yieldSize=1)))
  }
  sapply(bamfiles, getRL1)
}
bamCovToGRanges <- function(bamfiles, gr, lib.sizes, n=1000) {
  locs <- round(seq(from=1,to=width(gr),length=n))
  covmat <- sapply(seq_along(bamfiles), function(i) {
                     bf <- bamfiles[i]
                     cov <- coverage(bf, param=ScanBamParam(which=gr))[gr][[1]]
                     as.numeric(cov[locs]) * 1e6 / lib.sizes[i]
                   })
  gr.out <- GRanges(seqnames(gr), IRanges(locs + start(gr), width=1))
  mcols(gr.out) <- covmat
  gr.out
}
getGCOverGrid <- function(genome, gr, n=1000) {
  firstbp <- start(gr)
  dna <- getSeq(genome, gr)[[1]]
  j <- round(seq(from=1, to=length(dna)-1, length=n+1))
  start <- j[-(n+1)]
  end <- j[-1] + 1
  v <- Views(dna, start=start, end=end)
  gc <- as.numeric(letterFrequency(v, letters="GC", as.prob=TRUE))
  list(gc=gc, start=start, end=end)
}
densityToMean <- function(d) {
  delta <- d$x[2] - d$x[1]
  delta * sum(d$x * d$y)
}
alpineFlag <- function() scanBamFlag(isSecondaryAlignment=FALSE)
readGAlignAlpine <- function(bamfile, generange, manual=TRUE) {
  if (manual) {
    param <- ScanBamParam(which=generange, what=c("flag","mrnm","mpos"), flag=alpineFlag())
    gal <- readGAlignments(bamfile, use.names=TRUE, param=param)
    makeGAlignmentPairs(gal)
  } else {
    readGAlignmentPairs(bamfile,param=ScanBamParam(which=generange,flag=alpineFlag()))
  }
}
