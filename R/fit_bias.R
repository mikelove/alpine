#' Fit bias model over single-isoform genes
#'
#' This function estimates bias parameters for a single sample
#' over a set of single-isoform
#' genes. ~100 medium to highly expressed genes should be sufficient to
#' estimate the bias parameters robustly.
#' 
#' @param genes a GRangesList with the exons of different genes
#' @param bamfile a character string pointing to an indexed BAM file
#' @param fragtypes the output of \link{buildFragtypesFromExons}. must contain
#' the potential fragment types for the genes named in \code{genes}
#' @param genome a BSGenome object
#' @param models a list of character strings or formula describing the bias models, see vignette
#' @param readlength the read length
#' @param minsize the minimum fragment length to model
#' @param maxsize the maximum fragment length to model
#' @param zerotopos the ratio of zero counts (unobserved fragments) to positive counts
#' (observed fragments) that should be used in the Poisson GLM. The default value is 2.
#' The unobserved fragment types are then downsampled to match this ratio, and
#' up-weighted in the Poisson GLM.
#' @param speedglm logical, whether to use speedglm to estimate the coefficients.
#' Default is TRUE.
#'
#' @return
#'
#' @export
fitModelOverGenes <- function(genes, bamfile, fragtypes, genome,
                              models, readlength, minsize, maxsize,
                              zerotopos=2, speedglm=TRUE) {
  stopifnot(file.exists(bamfile))
  stopifnot(file.exists(paste0(as.character(bamfile),".bai")))
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
  fpbp <- getFPBP(genes, bamfile)
  # want ~1000 rows per gene, so ~300 reads per gene
  # so ~300/1500 = 0.2 fragments per basepair 
  target.fpbp <- 0.4
  fitpar.sub <- list()
  fitpar.sub[["coefs"]] <- list()
  # create a list over genes, populated with read info from this 'bamfile'
  # so we create a new object, and preserve the original 'fragtypes' object
  fragtypes.sub.list <- list()
  for (i in seq_along(genes)) {
    gene.name <- names(genes)[i]
    gene <- genes[[gene.name]]
    l <- sum(width(gene))
    # add counts per sample and subset
    generange <- range(gene)
    strand(generange) <- "*" # not necessary
    if (!as.character(seqnames(generange)) %in% seqlevels(BamFile(bamfile))) next
    # this necessary to avoid hanging on highly duplicated regions
    ## roughNumFrags <- countBam(bamfile, param=ScanBamParam(which=generange))$records/2
    ## if (roughNumFrags > 10000) next
    suppressWarnings({
                       ga <- readGAlignAlpine(bamfile, generange)
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
    # subset to include all (N) fragment locations and (N * zerotopos) zero locations
    # now we build a list of subsetted fragtypes
    fragtypes.sub.list[[gene.name]] <- subsetAndWeightFraglist(fraglist.temp, zerotopos)
  }
  if (length(fragtypes.sub.list) == 0) stop("not enough reads to model: ",bamfile)
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
    fivep <- fragtypes.sub$fivep[fragtypes.sub$fivep.test]
    threep <- fragtypes.sub$threep[fragtypes.sub$threep.test]
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
        fit <- glm(formula(f), family=poisson, data=fragtypes.sub, weights=wts, offset=offset)
      }
    fitpar.sub[["coefs"]][[modeltype]] <- fit$coefficients
    fitpar.sub[["summary"]][[modeltype]] <- summary(fit)$coefficients
  }
  fitpar.sub
}
