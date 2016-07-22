#' Fit bias models over single-isoform genes
#'
#' This function estimates parameters for one or more bias models
#' for a single sample over a set of single-isoform
#' genes. ~100 medium to highly expressed genes should be sufficient to
#' estimate the parameters robustly.
#' 
#' @param genes a GRangesList with the exons of different
#' single-isoform genes
#' @param bam.file a character string pointing to an indexed BAM file
#' @param fragtypes the output of \link{buildFragtypes}. must contain
#' the potential fragment types for the genes named in \code{genes}
#' @param genome a BSgenome object
#' @param models a list of lists: the outer list describes multiple models
#' each element of the inner list has two elements: \code{formula} and \code{offset}.
#' \code{formula} should be a character strings of an R formula
#' describing the bias models, e.g. \code{"count ~ ns(gc) + gene"}.
#' \code{offset} should be a character vector
#' listing possible bias offsets to be used (\code{"fraglen"} or \code{"vlmm"}).
#' Either \code{offset} or \code{formula} can be NULL for a model.
#' See vignette for recommendations and details.
#' @param readlength the read length
#' @param minsize the minimum fragment length to model
#' @param maxsize the maximum fragment length to model
#' @param speedglm logical, whether to use speedglm to estimate the coefficients.
#' Default is TRUE.
#'
#' @return a list with elements: coefs, summary, models,
#' and optional offets: fraglen.density, vlmm.fivep,
#' and vlmm.threep.
#' \itemize{
#' \item \strong{coefs} gives the estimated coefficients
#' for the different models that specified formula.
#' \item \strong{summary} gives the tables with coefficients, standard
#' errors and p-values,
#' \item \strong{models} stores the incoming
#' \code{models} list,
#' \item \strong{fraglen.density} is a
#' estimated density object for the fragment length distribution,
#' \item \strong{vlmm.fivep} and \strong{vlmm.threep}
#' store the observed and expected tabulations for the different
#' orders of the VLMM for read start bias.
#' }
#'
#' @references
#'
#' The complete bias model including fragment sequence bias:
#'
#' Love, M.I., Hogenesch, J.B., and Irizarry, R.A.,
#' Modeling of RNA-seq fragment sequence bias reduces
#' systematic errors in transcript abundance estimation.
#' bioRxiv (2015) doi: 10.1101/025767
#' 
#' The read start VLMM:
#'
#' Roberts, A., Trapnell, C., Donaghey, J., Rinn, J.L., and Pachter, L.,
#' Improving RNA-Seq expression estimates by correcting for fragment bias.
#' Genome Biology (2011) doi: 10.1186/gb-2011-12-3-r22
#' 
#' @examples
#'
#' # see vignette for a more realistic example
#' 
#' library(GenomicRanges)
#' library(BSgenome.Hsapiens.NCBI.GRCh38)
#' data(preprocessedData)
#' readlength <- 100
#' minsize <- 125 # see vignette how to choose
#' maxsize <- 175 # see vignette how to choose
#' # here a very small subset, should be ~100 genes
#' gene.names <- names(ebt.fit)[6:8] 
#' names(gene.names) <- gene.names
#' fragtypes <- lapply(gene.names, function(gene.name) {
#'                       buildFragtypes(ebt.fit[[gene.name]],
#'                                      Hsapiens, readlength,
#'                                      minsize, maxsize)
#' })
#' models <- list(
#'   "GC" = list(formula = "count ~ ns(gc,knots=gc.knots,
#'                 Boundary.knots=gc.bk) + gene",
#'               offset=c("fraglen"))
#' )
#' bam.file <- "~/proj/alpine/alpine/inst/extdata/ERR188088_galignpairs.bam"
#' fitpar <- fitBiasModels(genes=ebt.fit[gene.names],
#'                         bam.file=bam.file,
#'                         fragtypes=fragtypes,
#'                         genome=Hsapiens,
#'                         models=models,
#'                         readlength=readlength,
#'                         minsize=minsize,
#'                         maxsize=maxsize)
#'
#' @export
fitBiasModels <- function(genes, bam.file, fragtypes, genome,
                          models, readlength, minsize, maxsize,
                          speedglm=TRUE) {
  stopifnot(file.exists(bam.file))
  stopifnot(file.exists(paste0(as.character(bam.file),".bai")))
  stopifnot(is(genes, "GRangesList"))
  stopifnot(all(!is.na(sapply(models, function(x) x$formula))))
  stopifnot(is.numeric(readlength) & length(readlength) == 1)
  stopifnot(all(names(genes) %in% names(fragtypes)))
  if (any(sapply(models, function(m) "vlmm" %in% m$offset))) {
    stopifnot("fivep" %in% colnames(fragtypes[[1]]))
  }
  exon.dna <- getSeq(genome, genes)
  gene.seqs <- as(lapply(exon.dna, unlist), "DNAStringSet")
  # FPBP needed to downsample to a target fragment per kilobase
  fpbp <- getFPBP(genes, bam.file)
  # want ~1000 rows per gene, so ~300 reads per gene
  # so ~300/1500 = 0.2 fragments per basepair 
  target.fpbp <- 0.4
  fitpar.sub <- list()
  fitpar.sub[["coefs"]] <- list()
  fitpar.sub[["summary"]] <- list()
  # create a list over genes, populated with read info from this 'bam.file'
  # so we create a new object, and preserve the original 'fragtypes' object
  fragtypes.sub.list <- list()
  for (i in seq_along(genes)) {
    gene.name <- names(genes)[i]
    gene <- genes[[gene.name]]
    l <- sum(width(gene))
    # add counts per sample and subset
    generange <- range(gene)
    strand(generange) <- "*" # not necessary
    if (!as.character(seqnames(generange)) %in% seqlevels(BamFile(bam.file))) next
    # this necessary to avoid hanging on highly duplicated regions
    ## roughNumFrags <- countBam(bam.file, param=ScanBamParam(which=generange))$records/2
    ## if (roughNumFrags > 10000) next
    suppressWarnings({
                       ga <- readGAlignAlpine(bam.file, generange)
                     })
    if (length(ga) < 20) next
    ga <- keepSeqlevels(ga, as.character(seqnames(gene)[1]))
    # TODO: do this *after* finding compatible overlaps?
    # downsample to a target FPBP
    nfrags <- length(ga)
    this.fpbp <- nfrags / l
    if (this.fpbp > target.fpbp) {
      ga <- ga[sample(nfrags, round(nfrags * target.fpbp / this.fpbp), FALSE)]
    } 
    fco <- findCompatibleOverlaps(ga, GRangesList(gene))
    # message("-- ",round(length(fco)/length(ga),2)," compatible overlaps")
    # as.numeric(table(as.character(strand(ga))[queryHits(fco)])) # strand balance
    reads <- gaToReadsOnTx(ga, GRangesList(gene), fco)
    # fraglist.temp is a list of length 1
    # ...(matchReadsToFraglist also works for multiple transcripts)
    # it will only last for a few lines...
    fraglist.temp <- matchReadsToFraglist(reads, fragtypes[gene.name])
    # remove first and last bp for fitting the bias terms
    not.first.or.last.bp <- !(fraglist.temp[[1]]$start == 1 | fraglist.temp[[1]]$end == l)
    fraglist.temp[[1]] <- fraglist.temp[[1]][not.first.or.last.bp,]
    if (sum(fraglist.temp[[1]]$count) < 20) next
    # randomly downsample and up-weight
    fragtypes.sub.list[[gene.name]] <- subsetAndWeightFraglist(fraglist.temp,
                                                               downsample=200,
                                                               minzero=700)
  }
  if (length(fragtypes.sub.list) == 0) stop("not enough reads to model: ",bam.file)
  # collapse the list over genes into a
  # single DataFrame with the subsetted and weighted
  # potential fragment types from all genes
  # message("num genes w/ suf. reads: ",length(fragtypes.sub.list))
  if (length(fragtypes.sub.list) < 2) stop("requires at least two genes to fit model")
  gene.nrows <- sapply(fragtypes.sub.list, nrow)
  # message("mean rows per gene: ", round(mean(gene.nrows)))
  # a DataFrame of the subsetted fragtypes
  fragtypes.sub <- do.call(rbind, fragtypes.sub.list)

  # check the FPBP after downsampling:
  ## gene.counts <- sapply(fragtypes.sub.list, function(x) sum(x$count))
  ## gene.lengths <- sum(width(genes))
  ## round(unname(gene.counts / gene.lengths[names(gene.counts)]), 2)

  fitpar.sub[["models"]] <- models
  
  if (any(sapply(models, function(m) "fraglen" %in% m$offset))) {
    ## -- fragment bias --
    pos.count <- fragtypes.sub$count > 0
    fraglens <- rep(fragtypes.sub$fraglen[pos.count], fragtypes.sub$count[pos.count])
    fraglen.density <- density(fraglens)
    fragtypes.sub$logdfraglen <- log(matchToDensity(fragtypes.sub$fraglen, fraglen.density))
    # with(fragtypes.sub, plot(fraglen, exp(logdfraglen), cex=.1))
    fitpar.sub[["fraglen.density"]] <- fraglen.density
  }

  if (any(sapply(models, function(m) "vlmm" %in% m$offset))) {
    ## -- random hexamer priming bias with VLMM --
    pos.count <- fragtypes.sub$count > 0
    fivep <- fragtypes.sub$fivep[fragtypes.sub$fivep.test & pos.count]
    threep <- fragtypes.sub$threep[fragtypes.sub$threep.test & pos.count]
    vlmm.fivep <- fitVLMM(fivep, gene.seqs)
    vlmm.threep <- fitVLMM(threep, gene.seqs)
    ## par(mfrow=c(2,1))
    ## plotOrder0(vlmm.fivep$order0)
    ## plotOrder0(vlmm.threep$order0)
    
    # now calculate log(bias) for each fragment based on the VLMM
    fragtypes.sub <- addVLMMBias(fragtypes.sub, vlmm.fivep, vlmm.threep)
    fitpar.sub[["vlmm.fivep"]] <- vlmm.fivep
    fitpar.sub[["vlmm.threep"]] <- vlmm.threep
  }

  # allow a gene-specific intercept (although mostly handled already with downsampling)
  fragtypes.sub$gene <- factor(rep(seq_along(gene.nrows), gene.nrows))
  for (modeltype in names(models)) {
    if (is.null(models[[modeltype]]$formula)) {
      next
    }
    gc.knots <- seq(from=.4, to=.6, length=3)
    gc.bk <- c(0,1)
    relpos.knots <- seq(from=.25, to=.75, length=3)
    relpos.bk <- c(0,1)
    # message("fitting model type: ",modeltype)
    f <- models[[modeltype]]$formula
    offset <- numeric(nrow(fragtypes.sub))
    if ("fraglen" %in% models[[modeltype]]$offset) {
      # message("-- fragment length correction")
      offset <- offset + fragtypes.sub$logdfraglen
    }
    if ("vlmm" %in% models[[modeltype]]$offset) {
      # message("-- VLMM fragment start/end correction")
      offset <- offset + fragtypes.sub$fivep.bias + fragtypes.sub$threep.bias
    }
    if (!all(is.finite(offset))) stop("offset needs to be finite")
    fragtypes.sub$offset <-  offset
    if ( speedglm ) {
      # mm.small <- sparse.model.matrix(f, data=fragtypes.sub)
      mm.small <- model.matrix(formula(f), data=fragtypes.sub)
      stopifnot(all(colSums(abs(mm.small)) > 0))
      fit <- speedglm.wfit(fragtypes.sub$count, mm.small,
                           family=poisson(), 
                           weights=fragtypes.sub$wts,
                           offset=fragtypes.sub$offset)
    } else {
      fit <- glm(formula(f),
                 family=poisson,
                 data=fragtypes.sub,
                 weights=fragtypes.sub$wts,
                 offset=fragtypes.sub$offset)
    }
    fitpar.sub[["coefs"]][[modeltype]] <- fit$coefficients
    fitpar.sub[["summary"]][[modeltype]] <- summary(fit)$coefficients
  }
  fitpar.sub
}
