#' Estimate bias-corrected transcript abundances (FPKM) 
#'
#' This function takes the fitted bias parameters from \link{fitBiasModels}
#' and uses this information to derive bias corrected estimates of
#' transcript abundance for a gene (with one or more isoforms)
#' across multiple samples.
#' 
#' @param transcripts a GRangesList of the exons for multiple isoforms of a gene.
#' For a single-isoform gene, just wrap the exons in \code{GRangesList()}
#' @param bam.files a named vector pointing to the indexed BAM files
#' @param fitpar the output of \link{fitBiasModels}
#' @param genome a BSGenome object
#' @param models a list of character strings or formula describing the bias models, see vignette
#' @param readlength the read length
#' @param minsize the minimum fragment length to model
#' @param maxsize the maximum fragment length to model
#' @param subset logical, whether to downsample the non-observed fragments. Default is TRUE
#' @param niter the number of EM iterations. Default is 100.
#' @param lib.sizes a named vector of library sizes to use in calculating the FPKM.
#' If NULL (the default) a value of 1e6 is used for all samples.
#' @param optim logical, whether to use numerical optimization instead of the EM.
#' Default is FALSE.
#' @param custom.features an optional function to add custom features
#' to the fragment types DataFrame. This function takes in a DataFrame
#' returned by \link{buildFragtypes} and returns a DataFrame
#' with additional columns added. Default is NULL, adding no custom features.
#'
#' @return a list of lists. For each sample, a list with elements:
#' theta, lambda and count.
#' \itemize{
#' \item \strong{theta} gives the FPKM estimates for the
#' isoforms in \code{transcripts}
#' \item \strong{lambda} gives the average bias term
#' for the isoforms
#' \item \strong{count} gives the number of fragments which are
#' compatible with any of the isoforms in \code{transcripts}
#' }
#'
#' @examples
#'
#' # see vignette for a more realistic example
#'
#' # these next lines just write out a BAM file from R
#' # typically you would already have a BAM file
#' library(alpineData)
#' library(GenomicAlignments)
#' library(rtracklayer)
#' gap <- ERR188088()
#' dir <- system.file(package="alpineData", "extdata")
#' bam.file <- c("ERR188088" = file.path(dir,"ERR188088.bam"))
#' export(gap, con=bam.file)
#' 
#' data(preprocessedData)
#' library(GenomicRanges)
#' library(BSgenome.Hsapiens.NCBI.GRCh38)
#' models <- list(
#'  "GC"=list(formula="count~
#'  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
#'  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
#'  0",
#'  offset=c("fraglen"))
#' )
#'
#' readlength <- 75
#' minsize <- 125 # see vignette how to choose
#' maxsize <- 175 # see vignette how to choose
#' txs <- txdf.theta$tx_id[txdf.theta$gene_id == "ENSG00000198918"]
#' 
#' res <- estimateTheta(transcripts=ebt.theta[txs],
#'                      bam.files=bam.file,
#'                      fitpar=fitpar.small,
#'                      genome=Hsapiens,
#'                      models=models,
#'                      readlength=readlength,
#'                      minsize=minsize,
#'                      maxsize=maxsize)
#' 
#' @export
estimateTheta <- function(transcripts, bam.files, fitpar, genome,
                          models, readlength, minsize, maxsize,
                          subset=TRUE, niter=100, 
                          lib.sizes=NULL, optim=FALSE,
                          custom.features=NULL) {
  
  stopifnot(is(transcripts, "GRangesList"))
  stopifnot(length(transcripts) >= 1)
  singleiso <- length(transcripts) == 1
  stopifnot(!is.null(names(transcripts)))
  stopifnot(all(c("exon_rank","exon_id") %in% names(mcols(transcripts[[1]]))))

  stopifnot(!is.null(fitpar))
  stopifnot(all(!is.null(names(bam.files))))
  stopifnot(all(names(bam.files) %in% names(fitpar)))
  stopifnot(all(file.exists(bam.files)))

  stopifnot(all(file.exists(paste0(bam.files, ".bai"))))
  if (!is.null(lib.sizes)) stopifnot(all(names(bam.files) %in% names(lib.sizes)))
  if (is.null(lib.sizes)) {
    lib.sizes <- rep(1e6, length(bam.files))
    names(lib.sizes) <- names(bam.files)
  }
  
  w <- sum(width(transcripts))

  # is VLMM one of the offsets for any model
  any.vlmm <- any(sapply(models, function(m) "vlmm" %in% m$offset))

  # TODO: give better output for genes with smaller length than minsize
  if (min(w) <= minsize + 1) return(NULL)

  if (min(w) <= maxsize) {
    maxsize <- min(w)
  }
  
  # TODO: come up with a check on whether models is compatible with fitpar
  
  # this is a list of fragment types for each transcript
  st <- system.time({
    # TODO: could also save time by only doing GC stretches if necessary
    fraglist <- lapply(seq_along(transcripts), function(i) {
      out <- buildFragtypes(transcripts[[i]], genome, readlength,
                            minsize, maxsize, vlmm=any.vlmm)
      # optionally add more features to the fragment types DataFrame
      if (!is.null(custom.features)) {
        out <- custom.features(out)
      }
      out$tx <- names(transcripts)[i]
      out
    })
  })

  #message("building fragment types: ",round(unname(st[3]),1)," seconds")
  names(fraglist) <- names(transcripts)

  res <- lapply(seq_along(bam.files), function(i) {
    bam.file <- bam.files[i]
    bamname <- names(bam.file)
    txrange <- unlist(range(transcripts))
    strand(txrange) <- "*"
    generange <- range(txrange)
    
    #message("align reads to txs")
    suppressWarnings({
      ga <- readGAlignAlpine(bam.file, generange)
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
      st <- system.time({
        fragtypes <- subsetAndWeightFraglist(fraglist.temp)
      })
     #message("subset and weight fragment types: ", round(unname(st[3]),1), " seconds")
    } else {
        fragtypes <- do.call(rbind, fraglist.temp)
        # this is done in subsetAndWeightFraglist()
        fragtypes$genomic.id <- paste0(fragtypes$gstart,"-",fragtypes$gread1end,"-",
                                       fragtypes$gread2start,"-",fragtypes$gend)
    }

    # message("fragment bias")
    ## -- fragment bias --
    fraglen.density <- fitpar[[bamname]][["fraglen.density"]]
    fragtypes$logdfraglen <- log(matchToDensity(fragtypes$fraglen, fraglen.density))

    if (any.vlmm) {
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
      # this gives list output for one BAM file
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

    # NOTE: duplicated weights and bias are not the same for each tx.
    # The bias will often be identical for read start bias,
    # and very close for fragment length and fragment GC content given long reads.
    # It will not be so similar for relative position bias.
    # Zhonghui Xu points out: why not do the extra bookkeeping and
    # have the proper lambda-hat_ij fill out the A matrix.
    fragtypes.sub <- fragtypes[!duplicated(fragtypes$genomic.id),,drop=FALSE]
    stopifnot(all(fragtypes.sub$genomic.id == colnames(mat)))

    #message("run EM for models: ",paste(names(models), collapse=", "))
    n.obs <- fragtypes.sub$count
    
    # run EM for different models
    # this gives list output for one BAM file
    res.sub <- lapply(model.names, function(modeltype) {
      log.lambda <- getLogLambda(fragtypes, models, modeltype, fitpar, bamname)
      ## pred0 <- as.numeric(exp(log.lambda))
      ## pred <- pred0/mean(pred0)*mean(fragtypes.sub$count)
      ## boxplot(pred ~ factor(cut(fragtypes.sub$count,c(-1:10 + .5,20,Inf))), main=modeltype, range=0)
      if (is.null(lib.sizes)) {
        N <- mean(n.obs)
      } else {
        # TODO: in addition to the interval of considered lengths L
        # account for the triangle of fragments not in the count matrix
        N <- lib.sizes[bamname] / (1e9 * (maxsize - minsize))
      }
      
      # transcript-specific bias
      lambda.mat <- mat
      for (tx in names(transcripts)) {
        tx.id <- fragtypes$genomic.id[fragtypes$tx == tx]
        tx.idx <- match(tx.id, colnames(mat))
        lambda.mat[tx, tx.idx] <- exp(log.lambda[fragtypes$tx == tx])
      }

      wts <- if (subset) { fragtypes.sub$wts } else { 1 }

      # A also includes the library size
      A <- N * lambda.mat 
      theta <- runEM(n.obs, A, wts, niter, optim)

      # the average lambda for each transcript is stored in results
      lambda <- if (subset) {
        lambda.mat %*% wts / mat %*% wts
      } else {
        rowSums(lambda.mat) / rowSums(mat)
      }

      lambda <- as.numeric(lambda)
      names(lambda) <- names(transcripts)
      list(theta=theta, lambda=lambda)
    })
    return(c(res.sub, count=numCompatible))
  })
  names(res) <- names(bam.files)
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
