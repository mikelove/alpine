#' Predict coverage for a single-isoform gene
#'
#' Predict coverage for a single-isoform gene given
#' fitted bias parameters in a set of models,
#' and compare to the observed fragment coverage.
#'
#' Note that if the range between \code{minsize} and \code{maxsize}
#' does not cover most of the fragment length distribution, the
#' predicted coverage will underestimate the observed coverage.
#'
#' @param gene a GRangesList with the exons of different genes
#' @param bam.files a character string pointing to indexed BAM files
#' @param fitpar the output of running \link{fitBiasModels}
#' @param genome a BSgenome object
#' @param model.names a character vector listing the models,
#' see same argument in \link{estimateAbundance}
#'
#' @return a list with elements frag.cov, the observed fragment coverage
#' from the \code{bam.files} and pred.cov, a list with the predicted
#' fragment coverage for each of the \code{models}. 
#'
#' @examples
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
#' library(BSgenome.Hsapiens.NCBI.GRCh38)
#' 
#' model.names <- c("fraglen","fraglen.vlmm","GC","all")
#' 
#' pred.cov <- predictCoverage(gene=ebt.fit[["ENST00000379660"]],
#'                             bam.files=bam.file,
#'                             fitpar=fitpar.small,
#'                             genome=Hsapiens,
#'                             model.names=model.names)
#' 
#' # plot the coverage:
#' # note that, because [125,175] bp range specified in fitpar.small
#' # does not cover the fragment width distribution, the predicted curves
#' # will underestimate the observed. we correct here post-hoc
#' 
#' frag.cov <- pred.cov[["ERR188088"]][["frag.cov"]]
#' plot(frag.cov, type="l", lwd=3, ylim=c(0,max(frag.cov)*1.5))
#' for (i in seq_along(model.names)) {
#'   m <- model.names[i]
#'   pred <- pred.cov[["ERR188088"]][["pred.cov"]][[m]]
#'   lines(pred/mean(pred)*mean(frag.cov), col=i+1, lwd=3)
#' }
#' legend("topright", legend=c("observed",model.names),
#'        col=seq_len(length(model.names)+1), lwd=3)
#' 
#' @export
predictCoverage <- function(gene, bam.files, fitpar, genome, model.names) {
  stopifnot(is(gene, "GRanges"))
  stopifnot(!is.null(fitpar))
  stopifnot(all(names(bam.files) %in% names(fitpar)))
  if (is.null(names(bam.files))) {
    names(bam.files) <- seq_along(bam.files)
  }

  # pull out some model parameters
  stopifnot(all(c("readlength","minsize","maxsize","maxsize") %in%
                names(fitpar[[1]][["model.params"]])))
  readlength <- fitpar[[1]][["model.params"]][["readlength"]]
  minsize <- fitpar[[1]][["model.params"]][["minsize"]]
  maxsize <- fitpar[[1]][["model.params"]][["maxsize"]]

  # take model names and fitpar models and make the
  # models suitable for bias calculation
  models <- namesToModels(model.names, fitpar)
  
  fragtypes <- buildFragtypes(gene, genome, readlength=readlength,
                              minsize=minsize, maxsize=maxsize)
  res <- list()
  for (bamname in names(bam.files)) {
    # add counts
    bam.file <- bam.files[bamname]
    generange <- range(gene)
    strand(generange) <- "*" # not necessary
    suppressWarnings({
      ga <- readGAlignAlpine(bam.file, generange)
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
    fragtypes.temp$logdfraglen <- log(matchToDensity(fragtypes.temp$fraglen,
                                                     fraglen.density))
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

    ir <- IRanges(fragtypes.temp$start, fragtypes.temp$end)
    res[[bamname]]$l <- l
    res[[bamname]]$frag.cov <- frag.cov
    res[[bamname]]$pred.cov <- list()
    for (modeltype in names(models)) {
      # message("predicting model type: ",modeltype)
      log.lambda <- getLogLambda(fragtypes.temp, models, modeltype, fitpar, bamname)
      pred0 <- exp(log.lambda)
      pred <- pred0/mean(pred0)*mean(fragtypes.temp$count)
      res[[bamname]][["pred.cov"]][[modeltype]] <- coverage(ir, weight=pred)
    }
  }
  res
}
