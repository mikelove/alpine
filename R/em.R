#' Estimate bias-corrected transcript abundances (FPKM) 
#'
#' This function takes the fitted bias parameters from \link{fitBiasModels}
#' and uses this information to derive bias corrected estimates of
#' transcript abundance for a gene (with one or more isoforms)
#' across multiple samples.
#' 
#' @param transcripts a GRangesList of the exons for multiple isoforms of a gene.
#' For a single-isoform gene, just wrap the exons in \code{GRangesList()}
#' @param bamfiles a named vector pointing to the indexed BAM files
#' @param fitpar the output of \link{fitBiasModels}
#' @param genome a BSGenome object
#' @param models a list of character strings or formula describing the bias models, see vignette
#' @param readlength the read length
#' @param minsize the minimum fragment length to model
#' @param maxsize the maximum fragment length to model
#' @param subset logical, whether to downsample the non-observed fragments. Default is TRUE
#' @param zerotopos the rate of downsampling, see \link{fitBiasModels}.
#' Here it is recommended to use a higher value than for fitting the bias parameters.
#' Default is 20.
#' @param niter the number of EM iterations. Default is 100.
#' @param lib.sizes a named vector of library sizes to use in calculating the FPKM.
#' If NULL (the default) a value of average coverage will be calculated from
#' the observed fragments. Better is to provide a named vector with 1e6 for all samples.
#' @param optim logical, whether to use numerical optimization instead of the EM.
#' Default is FALSE.
#'
#' @return a list of lists. For each sample, a list with elements:
#' theta, lambda and count. theta gives the FPKM estimates for the
#' isoforms in \code{transcripts}. lambda gives the average bias term
#' for the isoforms, and count gives the number of fragments which are
#' compatible with any of the isoforms in \code{transcripts}
#'
#' @export
estimateTheta <- function(transcripts, bamfiles, fitpar, genome,
                          models, readlength, minsize, maxsize,
                          subset=TRUE, zerotopos=20, niter=100, 
                          lib.sizes=NULL, optim=FALSE) {
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
      out <- buildFragtypes(transcripts[[i]], genome, readlength, minsize, maxsize)
      out$tx <- names(transcripts)[i]
      out
    })
  })
  # message("building fragment types: ",round(unname(st[3]))," seconds")
  names(fraglist) <- names(transcripts)

  res <- lapply(seq_along(bamfiles), function(i) {
    bamfile <- bamfiles[i]
    bamname <- names(bamfile)
    txrange <- unlist(range(transcripts))
    strand(txrange) <- "*"
    generange <- range(txrange)
    
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
      # this is also done in subsetAndWeightFraglist()
      fragtypes$genomic.id <- paste0(fragtypes$gstart,"-",fragtypes$gread1end,"-",
                                     fragtypes$gread2start,"-",fragtypes$gend)
    }

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

    fragtypes.sub <- fragtypes[!duplicated(fragtypes$genomic.id),,drop=FALSE]
    stopifnot(all(fragtypes.sub$genomic.id == colnames(mat)))

    #message("run EM for models: ",paste(names(models), collapse=", "))
    n.obs <- fragtypes.sub$count
    
    # run EM for different models
    # this gives list output for one bamfile
    res.sub <- lapply(model.names, function(modeltype) {

      # this vector goes over 'fragtypes' rows, so includes duplicates
      log.lambda <- getLogLambda(fragtypes, models, modeltype, fitpar, bamname)

      ## pred0 <- as.numeric(exp(log.lambda))
      ## pred <- pred0/mean(pred0)*mean(fragtypes.sub$count)
      ## boxplot(pred ~ factor(cut(fragtypes.sub$count,c(-1:10 + .5,20,Inf))),
      ##         main=modeltype, range=0)
      N <- if (is.null(lib.sizes)) {
        mean(n.obs)
      } else {
          # TODO: in addition to the interval of considered lengths L
          # account for the triangle of fragments not in the count matrix
        lib.sizes[bamname] / (1e9 * (maxsize - minsize))
      }
      
      # transcript-specific bias terms
      A <- mat * N
      for (tx in names(transcripts)) {
        tx.id <- fragtypes$genomic.id[fragtypes$tx == tx]
        tx.idx <- match(tx.id, colnames(A))
        A[tx, tx.idx] <- A[tx, tx.idx] * exp(log.lambda[fragtypes$tx == tx])
      }
      wts <- if (subset) { fragtypes.sub$wts } else { 1 }
      theta <- runEM(n.obs, A, wts, niter, optim)

      # TODO fix this....
      
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

######### unexported EM functions #########

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
