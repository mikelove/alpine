incidenceMat <- function(x, y, numeric=TRUE) {
  # borrowed from Wolfgang Huber
  ux = unique(x)
  uy = unique(y)
  im = matrix(FALSE, nrow=length(ux), ncol=length(uy), dimnames=list(ux, uy))
  im[ cbind(x, y) ] = TRUE
  if (numeric) {
    mode(im) <- "numeric"
  }
  return(im)
}
runEM <- function(n.obs, A, wts=1, niter=20, optim=FALSE) {
  J <- ncol(A)
  ntx <- nrow(A)
  log.like <- function(theta.hat) {
    sum(wts * dpois(n.obs, colSums(A * theta.hat), log=TRUE))
  }
  theta.hat <- rep(1, ntx)
  theta.0 <- rep(1, ntx)
  n.obs.sub <- n.obs[n.obs > 0]
  A.sub <- A[,n.obs > 0,drop=FALSE]
  rowSumsA <- rowSums(t(t(A) * wts))
  if (!optim) {
    for (tt in 1:niter) {
      n.hat <- t(t(theta.hat * A.sub) * n.obs.sub / colSums(theta.hat * A.sub))
      theta.hat <- rowSums(n.hat) / rowSumsA
    }
  } else {
    theta.hat <- optim(theta.hat, log.like,
                       lower=rep(1e-6,ntx), upper=rep(1e6,ntx),
                       control=list(fnscale=-1), method="L-BFGS-B")$par
  }
  theta.hat
}
estimateTheta <- function(transcripts, bamfiles, fitpar, genome,
                          models, readlength,
                          subset=FALSE, zerotopos, niter, 
                          lib.sizes=NULL, optim=FALSE, minsize, maxsize) {
  stopifnot(is(transcripts, "GRangesList"))
  stopifnot(length(transcripts) >= 1)
  singleiso <- length(transcripts) == 1
  stopifnot(all(c("exon_rank","exon_id") %in% names(mcols(transcripts[[1]]))))
  stopifnot(all(names(bamfiles) %in% names(fitpar)))
  stopifnot(all(!is.null(names(bamfiles))))
  stopifnot(all(file.exists(bamfiles)))
  stopifnot(!is.null(fitpar))
  stopifnot(all(names(bamfiles) %in% names(fitpar)))
  stopifnot(!is.null(names(transcripts)))
  stopifnot(all(file.exists(paste0(bamfiles, ".bai"))))
  if (!is.null(lib.sizes)) stopifnot(all(names(bamfiles) %in% names(lib.sizes)))
  w <- sum(width(transcripts))

  # TODO: give better output here
  if (min(w) <= minsize + 1) return(NULL)

  if (min(w) <= maxsize) {
    maxsize <- min(w)
  }
  
  # TODO: come up with a check on whether models is compatible with fitpar
  
  # this is a list of fragment types for each transcript
  st <- system.time({ 
    fraglist <- lapply(seq_along(transcripts), function(i) {
      out <- buildFragtypesFromExons(transcripts[[i]], genome, readlength, minsize, maxsize)
      out$tx <- names(transcripts)[i]
      out
    })
  })
  # message("building fragment types: ",round(unname(st[3]))," seconds")
  names(fraglist) <- names(transcripts)
  res <- lapply(seq_along(bamfiles), function(i) {
    bamfile <- bamfiles[i]
    bamname <- names(bamfile)
    generange <- range(unlist(range(transcripts)))
    strand(generange) <- "*" # not necessary
    #message("align reads to txs")
    suppressWarnings({
      ga <- readGAlignAlpine(bamfile, generange)
    })
    #message("-- ",length(ga)," reads")

    outputZero <- FALSE
    if (length(ga) == 0) {
      numCompatible <- 0
      outputZero <- TRUE
    } else { 
      ga <- keepSeqlevels(ga, as.character(seqnames(transcripts[[1]])[1]))
      fco <- findCompatibleOverlaps(ga, transcripts)
      numCompatible <- length(unique(queryHits(fco)))
      #message("-- ",round(numCompatible/length(ga),2)," compatible overlaps")
      #message("---- ",seqnames(generange),":",start(generange),"-",end(generange))
      # table(strand(ga)[unique(queryHits(fco))]) # are the read counts even across strand?
      # boxplot(lapply(reads, function(x) width(x)))
      # here called "reads" although they are fragments. everything is already called fragments :-/
      reads <- gaToReadsOnTx(ga, transcripts, fco)
      fraglist.temp <- matchReadsToFraglist(reads, fraglist)
      txcounts <- sapply(fraglist.temp, function(x) sum(x$count))
      #message("---- ",paste(txcounts, collapse=" "))
      if (all(txcounts == 0)) outputZero <- TRUE
    }

    # report 0 output for all models if all txs have 0 count
    model.names <- names(models)
    names(model.names) <- model.names
    nms.tx <- names(transcripts)
    if (outputZero) {
      #message("all transcripts have 0 counts")
      theta <- numeric(length(nms.tx))
      lambda <- rep(NA,length(nms.tx)) # don't bother calculating lambda
      names(lambda) <- names(theta) <- nms.tx
      # for all models:
      res.sub <- lapply(model.names, function(x) {
                   list(theta=theta, lambda=lambda)
                 })
      return(c(res.sub,count=0)) # return results for this sample
    }
    
    if (subset) {
      #message("subset and weight fragment types")
      fragtypes <- subsetAndWeightFraglist(fraglist.temp, zerotopos)
    } else {
      fragtypes <- do.call(rbind, fraglist.temp)
    }
    fragtypes$genomic.id <- paste0(fragtypes$gstart,"-",fragtypes$gread1end,"-",
                                   fragtypes$gread2start,"-",fragtypes$gend)

    # message("fragment bias")
    ## -- fragment bias --
    fraglen.density <- fitpar[[bamname]][["fraglen.density"]]
    fragtypes$logdfraglen <- log(matchToDensity(fragtypes$fraglen, fraglen.density))

    if (any(sapply(models, function(m) "vlmm" %in% m$offset))) {
      stopifnot( "vlmm.fivep" %in% names(fitpar[[bamname]]) )
      # message("priming bias")
      ## -- random hexamer priming bias with VLMM --
      vlmm.fivep <- fitpar[[bamname]][["vlmm.fivep"]]
      vlmm.threep <- fitpar[[bamname]][["vlmm.threep"]]
      fragtypes <- addVLMMBias(fragtypes, vlmm.fivep, vlmm.threep)
    }
      
    # specific code for one isoform
    if (singleiso) {
      n.obs <- fragtypes$count
      # this gives list output for one bamfile
      res.sub <- lapply(model.names, function(modeltype) {
        log.lambda <- getLogLambda(fragtypes, models, modeltype, fitpar, bamname)
        log.lambda <- as.numeric(log.lambda)
        N <- if (is.null(lib.sizes)) {
          mean(n.obs)
        } else {
            # TODO: here fix like below
            lib.sizes[bamname] / (1e9 * (maxsize - minsize))
          }
        A <- N * exp(log.lambda)
        wts <- if (subset) { fragtypes$wts } else { 1 } 
        theta <- sum(n.obs  * wts)/sum(A * wts)
        lambda <- sum(wts * exp(log.lambda)) / sum(wts)
        names(lambda) <- names(theta) <- names(transcripts)
        list(theta=theta, lambda=lambda)
      })
      return(c(res.sub,count=numCompatible)) # return results for this sample
    }

    # make incidence matrix
    # duplicate genomic ID across tx will be a single column 
    mat <- incidenceMat(fragtypes$tx, fragtypes$genomic.id)
    # make sure the rows are in correct order
    stopifnot(all(rownames(mat) == names(transcripts)))
    # NOTE: duplicated weights are not the same for each tx
    fragtypes.sub <- fragtypes[!duplicated(fragtypes$genomic.id),,drop=FALSE]
    stopifnot(all(fragtypes.sub$genomic.id == colnames(mat)))
    #message("run EM for models: ",paste(names(models), collapse=", "))
    n.obs <- fragtypes.sub$count
    
    # run EM for different models
    # this gives list output for one bamfile
    res.sub <- lapply(model.names, function(modeltype) {
      log.lambda <- getLogLambda(fragtypes.sub, models, modeltype, fitpar, bamname)
      log.lambda <- as.numeric(log.lambda)
      ## pred0 <- as.numeric(exp(log.lambda))
      ## pred <- pred0/mean(pred0)*mean(fragtypes.sub$count)
      ## boxplot(pred ~ factor(cut(fragtypes.sub$count,c(-1:10 + .5,20,Inf))), main=modeltype, range=0)
      N <- if (is.null(lib.sizes)) {
        mean(n.obs)
      } else {
          # TODO: in addition to the interval of considered lengths L
          # account for the triangle of fragments not in the count matrix
        lib.sizes[bamname] / (1e9 * (maxsize - minsize))
      }
      A <- t(t(mat * N) * exp(log.lambda))
      wts <- if (subset) { fragtypes.sub$wts } else { 1 }
      theta <- runEM(n.obs, A, wts, niter, optim)
      lambda <- mat %*% (wts * exp(log.lambda)) / mat %*% wts
      lambda <- as.numeric(lambda)
      names(lambda) <- names(transcripts)
      list(theta=theta, lambda=lambda)
    })
    return(c(res.sub, count=numCompatible))
  })
  names(res) <- names(bamfiles)
  res
}
