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
# plot for fragment GC bias over multiple samples
# This takes the natural spline estimated coefficients from the Poisson regression
# and generates a smooth function of log fragment rate over fragment GC.
# This is similar to the plot you would get with, e.g. plot.gam()
# after having fit a generalized additive model.
plotGC <- function(fitpar, knots=c(.4,.5,.6), bk=c(0,1), col, lty, m="GC") {
  n <- length(knots)
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[m]][1:(n+2)])
  genecoefs <- lapply(fitpar, function(elem) elem[["coefs"]][[m]][ grep("gene", names(elem[["coefs"]][[m]])) ])
  # new intercept: the average of the intercept + gene coefficients
  coefmat[1,] <- coefmat[1,] + sapply(genecoefs, mean)
  z <- seq(from=.2,to=.8,length=101)
  x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk))
  logpred <- x %*% coefmat
  plot(0,0,type="n",xlim=c(.2,.8),ylim=c(min(logpred),max(logpred)),
       ylab="log fragment rate", xlab="fragment GC content",
       main="fragment sequence effect")
  if (missing(col)) {
    col <- rep("black", ncol(logpred))
  }
  if (missing(lty)) {
    lty <- rep(1, ncol(logpred))
  }
  for (i in 1:ncol(logpred)) {
    lines(z, logpred[,i], col=col[i], lwd=2, lty=lty[i])
  }
}
# plots for VLMM (read start sequence bias) for a single sample
plotOrder0 <- function(order0, ...) {
  dna.letters <- c("A","C","G","T")
  mat <- log(order0$obs/order0$expect)
  xpos <- -8:12
  dna.cols <- c("green3","blue3","orange3","red3")
  plot(0,0,xlim=c(-8,12),ylim=c(-0.3,0.3),type="n",xlab="position",ylab="log(observed / expected)",...)
  for (i in 1:4) {
    points(xpos, mat[i,], col=dna.cols[i], type="b", lwd=2)
  }
  abline(v=0, h=0, col=rgb(0,0,0,.3))
  legend("topright",dna.letters,pch=1,lty=1,col=dna.cols,bg="white")
}
plotOrder1 <- function(order1, pos1) {
  order <- 1
  npos <- length(pos1)
  dna.cols <- c("green3","blue3","orange3","red3")
  mypar(1,npos+1,mar=c(5,1,3,1))
  for (i in 1:npos) {
    plot(as.vector(log(order1$obs[,,i]/order1$expect)), rev(seq_len(4 * 4^order)),
         col=rep(dna.cols,each=4), xlim=c(-1,1),
         ylab="",xlab="",yaxt="n",main=pos1[i] - 10 + 1,cex=2)
    abline(v=0)
  }
  plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  alpha <- alphafun(dna.letters, order-1)
  legend("center",alpha,pch=1,col=dna.cols,cex=2,title="prev")
}
plotOrder2 <- function(order2, pos2) {
  order <- 2
  npos <- length(pos2)
  dna.cols <- c("green3","blue3","orange3","red3")
  mypar(1,npos+1,mar=c(5,1,3,1))
  for (i in 1:npos) {
    plot(as.vector(log(order2$obs[,,i]/order2$expect)), rev(seq_len(4 * 4^order)),
         col=rep(dna.cols,each=4),pch=rep(1:4,each=16), xlim=c(-1,1),
         ylab="",xlab="",yaxt="n",main=pos2[i] - 10 + 1,cex=2)
    abline(v=0)
  }
  plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  alpha <- alphafun(dna.letters, order-1)
  legend("center",alpha,pch=rep(1:4,each=4),col=dna.cols,cex=2,title="prev")
}
plotRelPos <- function(fitpar, knots=c(.25,.5,.75), bk=c(0,1), col, lty, m="GC") {
  n <- length(knots)
  # assuming same number of knots for GC and for relpos
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[m]][c(1,(3+n):(3+2*n))])
  z <- seq(from=0,to=1,length=101)
  x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk))
  logpred <- x %*% coefmat
  logpred <- scale(logpred, scale=FALSE)
  plot(0,0,type="n",xlim=c(0,1),ylim=c(min(logpred),max(logpred)),
       ylab="log fragment rate", xlab="5' -- position in transcript -- 3'",
       main="relative position bias")
  if (missing(col)) {
    col <- rep("black", ncol(logpred))
  }
  if (missing(lty)) {
    lty <- rep(1, ncol(logpred))
  }
  for (i in 1:ncol(logpred)) {
    lines(z, logpred[,i], col=col[i], lwd=2, lty=lty[i])
  }
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
readGAlignAlpine <- function(bamfile, generange) {
  readGAlignmentPairs(bamfile,param=ScanBamParam(which=generange,flag=alpineFlag()))
}
