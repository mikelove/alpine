#' alpine: bias corrected transcript abundance estimation
#'
#' alpine is a package for estimating and visualizing many forms of sample-specific
#' biases that can arise in RNA-seq, including fragment length
#' distribution, positional bias on the transcript, read
#' start bias (random hexamer priming), and fragment GC content
#' (amplification). It also offers bias-corrected estimates of
#' transcript abundance (FPKM). It is currently designed for
#' un-stranded paired-end RNA-seq data, but will be extended in the
#' future to support strand-specific data.
#'
#' See the vignette for a detailed workflow.
#' 
#' The main functions in this package are:
#' \enumerate{
#' \item \link{buildFragtypes} - build out features for fragment types from exons of a single gene (GRanges)
#' \item \link{fitBiasModels} - fit parameters for one or more bias models over a set of ~100 medium to highly expressed single isoform genes (GRangesList)
#' \item \link{estimateTheta} - given a set of genome alignments (BAM files) and a set of isoforms of a gene (GRangesList), estimate the transcript abundances for these isoforms (FPKM) for various bias models
#' \item \link{extract} - given a list of output from \code{estimateTheta}, compile an FPKM matrix across transcripts and samples
#' }
#'
#' Some helper functions for preparing gene objects:
#' \enumerate{
#' \item \link{splitGenesAcrossChroms} - split apart "genes" where isoforms are on different chromosomes
#' \item \link{splitLongGenes} - split apart "genes" which cover a suspiciously large range, e.g. 1 Mb
#' \item \link{mergeGenes} - merge overlapping isoforms into new "genes"
#' }
#' 
#' The plotting functions are:
#' \enumerate{
#' \item \link{plotGC} - plot the fragment GC bias curves
#' \item \link{plotFragLen} - plot the framgent length distributions
#' \item \link{plotRelPos} - plot the positional bias (5' to 3')
#' \item \link{plotOrder0}, \link{plotOrder1}, \link{plotOrder2} - plot the read start bias
#' }
#'
#' @references
#'
#' Michael I Love, John B Hogenesch, Rafael A Irizarry:
#' Modeling of RNA-seq fragment sequence bias reduces
#' systematic errors in transcript abundance estimation
#' Posted to bioRxiv August 2015, \url{http://biorxiv.org/content/early/2015/08/28/025767}
#'
#' @author Michael Love
#'
#' @importFrom IRanges match %in%
#' @importFrom splines ns
#' @importFrom speedglm speedglm.wfit
#' 
#' @docType package
#' @name alpine-package
#' @aliases alpine-package
#' @keywords package
NULL

#' Build fragment types from exons
#'
#' This is a core function used to construct a table of features used for
#' bias modeling, with one row for every potential fragment that could
#' arise from a transcript. The output of this function is used by
#' \link{fitBiasModels}, and this function is used inside \link{estimateTheta}
#' in order to model the bias affecting different fragments across isoforms
#' of a gene.
#' 
#' @param exons a GRanges object with the exons for a single transcript
#' @param genome a BSgenome object, for example \code{Hsapiens}, after
#' loading the package \code{BSgenome.Hsapiens.UCSC.hg19}
#' @param readlength the length of the reads. This doesn't necessarily
#' have to be exact (+/- 1 bp would be fine), but it should be close.
#' @param minsize the minimum fragment length to model. The interval between
#' \code{minsize} and \code{maxsize} should contain the central 95-99 percent
#' of the fragment length distribution
#' @param maxsize the maximum fragment length to model
#' @param gc logical, whether to estimate the fragment GC content
#' @param gc.str logical, whether to look for stretches of very high GC within fragments
#' @param vlmm logical, whether to estimate the Cufflinks variable length
#' Markov model for read starts
#'
#' @return a DataFrame with bias features for all potential fragments
#' 
#' @export
buildFragtypes <- function(exons, genome, readlength,
                                    minsize, maxsize, 
                                    gc=TRUE, gc.str=TRUE, vlmm=TRUE) {
  stopifnot(is(exons,"GRanges"))
  stopifnot(is(genome,"BSgenome"))
  stopifnot(is.numeric(minsize) & is.numeric(maxsize) & is.numeric(readlength))
  stopifnot(sum(width(exons)) >= maxsize)
  stopifnot(all(c("exon_rank","exon_id") %in% names(mcols(exons))))

  # these parameters must be fixed, as dictated by fitVLMM()
  npre <- 8
  npost <- 12
  
  map <- mapTxToGenome(exons)
  l <- nrow(map)
  strand <- as.character(strand(exons)[1])
  start <- rep(seq_len(l-minsize+1),each=maxsize-minsize+1)
  end <- as.integer(start + minsize:maxsize - 1)
  mid <- as.integer(0.5 * (start + end))
  relpos <- mid/l
  fraglen <- as.integer(end - start + 1)
  id <- IRanges(start, end)
  fragtypes <- DataFrame(start=start,end=end,relpos=relpos,fraglen=fraglen,id=id)
  fragtypes <- fragtypes[fragtypes$end <= l,,drop=FALSE]
  exon.dna <- getSeq(genome, exons)
  tx.dna <- unlist(exon.dna)
  if (vlmm) {
    # strings needed for VLMM
    fragtypes$fivep.test <- fragtypes$start - npre >= 1
    fragtypes$fivep <- as(Views(tx.dna, fragtypes$start - ifelse(fragtypes$fivep.test, npre, 0),
                                fragtypes$start + npost), "DNAStringSet")
    fragtypes$threep.test <- fragtypes$end + npre <= length(tx.dna) 
    fragtypes$threep <- as(Views(tx.dna, fragtypes$end - npost,
                                 fragtypes$end + ifelse(fragtypes$threep.test, npre, 0),),
                           "DNAStringSet")
    # reverse complement the three prime sequence
    fragtypes$threep <- reverseComplement(fragtypes$threep)
  }
  if (gc) {
    # get the GC content for the entire fragment
    fragrange <- minsize:maxsize
    gc.vecs <- lapply(fragrange, function(i) {
                        letterFrequencyInSlidingView(tx.dna, view.width=i, letters="CG", as.prob=TRUE)
                      })
    fragtypes <- fragtypes[order(fragtypes$fraglen),,drop=FALSE]
    fragtypes$gc <- do.call(c, gc.vecs)
    fragtypes <- fragtypes[order(fragtypes$start),,drop=FALSE]
  }
  if (gc.str) {
    # additional features: GC in smaller sections
    gc.40 <- as.numeric(letterFrequencyInSlidingView(tx.dna, 40, letters="CG", as.prob=TRUE))
    max.gc.40 <- max(Views(gc.40, fragtypes$start, fragtypes$end - 40 + 1))
    gc.20 <- as.numeric(letterFrequencyInSlidingView(tx.dna, 20, letters="CG", as.prob=TRUE))
    max.gc.20 <- max(Views(gc.20, fragtypes$start, fragtypes$end - 20 + 1))
    fragtypes$GC40.90 <- as.numeric(max.gc.40 >= 36/40)
    fragtypes$GC40.80 <- as.numeric(max.gc.40 >= 32/40)
    fragtypes$GC20.90 <- as.numeric(max.gc.20 >= 18/20)
    fragtypes$GC20.80 <- as.numeric(max.gc.20 >= 16/20)
  }
  # these are the fragment start and end in genomic space
  # so for minus strand tx, gstart > gend
  fragtypes$gstart <- txToGenome(fragtypes$start, map)
  fragtypes$gend <- txToGenome(fragtypes$end, map)
  fragtypes$gread1end <- txToGenome(fragtypes$start + readlength - 1, map)
  fragtypes$gread2start <- txToGenome(fragtypes$end - readlength + 1, map)
  #message("nrow fragtypes: ",nrow(fragtypes))
  fragtypes
}

######### unexported core functions #########

startLeft <- function(x) {
  first.plus <- as.logical(strand(first(x)) == "+")
  ifelse(first.plus, start(first(x)), start(last(x)))
}
endRight <- function(x) {
  first.plus <- as.logical(strand(first(x)) == "+")
  ifelse(first.plus, end(last(x)), end(first(x)))
}
mapTxToGenome <- function(exons) {
  strand <- as.character(strand(exons)[1])
  stopifnot(all(exons$exon_rank == seq_along(exons)))
  bases <- if (strand == "+") {
    do.call(c,lapply(exons, function(exon) start(exon):end(exon)))
  } else if (strand == "-") {
    do.call(c,lapply(exons, function(exon) end(exon):start(exon)))    
  }
  data.frame(tx=seq_along(bases),
             genome=bases,
             exon_rank=rep(exons$exon_rank, width(exons)))
}
genomeToTx <- function(genome, map) map$tx[match(genome, map$genome)]
txToGenome <- function(tx, map) map$genome[match(tx, map$tx)]
txToExon <- function(tx, map) map$exon_rank[match(tx, map$tx)]
gaToReadsOnTx <- function(ga, grl, fco=NULL) {
  reads <- list()
  for (i in seq_along(grl)) {
    exons <- grl[[i]]
    strand <- as.character(strand(exons)[1])
    read.idx <- if (is.null(fco)) {
      seq_along(ga)
    } else {
      queryHits(fco)[subjectHits(fco) == i]
    }
    map <- mapTxToGenome(exons)
    # depending on strand of gene:
    # start of left will be the first coordinate on the transcript (+ gene)
    # or start of left will be the last coordinate on the transcript (- gene)
    if (strand == "+") {
      start <- genomeToTx(startLeft(ga[read.idx]), map)
      end <- genomeToTx(endRight(ga[read.idx]), map)
    } else if (strand == "-") {
      start <- genomeToTx(endRight(ga[read.idx]), map)
      end <- genomeToTx(startLeft(ga[read.idx]), map)
    }
    valid <- start < end & !is.na(start) & !is.na(end)
    reads[[i]] <- IRanges(start[valid], end[valid])
  }
  names(reads) <- names(grl)
  reads
}
matchReadsToFraglist <- function(reads, fraglist) {
  for (tx.idx in seq_along(fraglist)) {
    unique.reads <- unique(reads[[tx.idx]])
    readtab <- table(match(reads[[tx.idx]], unique.reads))
    fraglist[[tx.idx]]$count <- 0
    reads.in.fraglist <- unique.reads %in% fraglist[[tx.idx]]$id
    unique.reads <- unique.reads[reads.in.fraglist]
    readtab <- readtab[reads.in.fraglist]    
    fraglist[[tx.idx]][match(unique.reads, fraglist[[tx.idx]]$id),"count"] <- as.numeric(readtab)
  }
  fraglist
}

subsetAndWeightFraglist <- function(fraglist, zerotopos, minzero=2000, maxmult=20000) {
  ntx <- length(fraglist)
  fraglist.sub <- list()
  for (tx.idx in seq_len(ntx)) {
    count <- fraglist[[tx.idx]]$count
    sumpos <- sum(count > 0)
    sumzero <- sum(count == 0)
    ## how many zeros to sample?
    # a multiple 'zerotopos' of the number of positive counts
    # or a preset minimum value 'minzero'
    multzero <- max(round(zerotopos*sumpos), minzero)
    # if this multiple is more than 'maxmult', we have plenty of zero,
    # so then sample either 'maxmult' or the number of positive counts,
    # whichever is larger
    if (multzero > maxmalt) {
      multzero <- max(maxmult, sumpos)
    }
    # of course, we cannot sample more than the total number of zero
    # 'numzero' is then the number of zero we will sample
    numzero <- min(sumzero, multzero)
    # down-sample
    idx <- c(which(count > 0), sample(which(count == 0), numzero, replace=FALSE))
    fraglist.sub[[tx.idx]] <- fraglist[[tx.idx]][idx,,drop=FALSE]
  }
  fragtypes <- do.call(rbind, fraglist.sub)
  fragtypes$genomic.id <- paste0(fragtypes$gstart,"-",fragtypes$gread1end,"-",
                                 fragtypes$gread2start,"-",fragtypes$gend)

  ######################### TODO
  zero.wt <- 20
  
  # return fragtypes, but with duplicate rows for selected fragments
  fragtypes <- fragtypes[fragtypes$genomic.id %in% fragtypes.sub$genomic.id,,drop=FALSE]
  fragtypes$wts <- rep(1, nrow(fragtypes))
  fragtypes$wts[fragtypes$count == 0] <- zero.wt
  fragtypes
}


matchToDensity <- function(x, d) {
  idx <- cut(x, c(-Inf, d$x, Inf))
  pdf <- c(0, d$y)
  pdf.x <- pdf[ idx ] + 1e-6
  stopifnot(all(pdf.x > 0))
  pdf.x
}
getFPBP <- function(genes, bamfile) {
  gene.ranges <- unlist(range(genes))
  gene.lengths <- sum(width(genes))
  res <- countBam(bamfile, param=ScanBamParam(which=gene.ranges))
  # two records per fragment
  out <- (res$records / 2)/gene.lengths
  names(out) <- names(genes)
  out
}
getLogLambda <- function(fragtypes, models, modeltype, fitpar, bamname) {
  # knots and boundary knots should come from the environ where the formula were defined
  f <- models[[modeltype]]$formula
  offset <- numeric(nrow(fragtypes))
  if ("fraglen" %in% models[[modeltype]]$offset) {
    # message("-- fragment length correction")
    offset <- offset + fragtypes$logdfraglen
  }
  if ("vlmm" %in% models[[modeltype]]$offset) {
    # message("-- VLMM fragment start/end correction")
    offset <- offset + fragtypes$fivep.bias + fragtypes$threep.bias
  }
  if (!is.null(f)) {
    gc.knots <- seq(from=.4, to=.6, length=3)
    gc.bk <- c(0,1)
    relpos.knots <- seq(from=.25, to=.75, length=3)
    relpos.bk <- c(0,1)
    # assume: no intercept in formula
    # sparse.model.matrix produces different column names, so don't use
    # mm.big <- sparse.model.matrix(f, data=fragtypes)
    mm.big <- model.matrix(formula(f), data=fragtypes)
    beta <- fitpar[[bamname]][["coefs"]][[modeltype]]
    stopifnot(any(colnames(mm.big) %in% names(beta)))
    if (all(is.na(beta))) stop("all coefs are NA")
    beta[is.na(beta)] <- 0 # replace NA coefs with 0: these were not observed in the training data
    # this gets rid of the gene1, gene2 and Intercept terms
    beta <- beta[match(colnames(mm.big), names(beta))]
    # add offset
    log.lambda <- as.numeric(mm.big %*% beta) + offset
  } else {
    log.lambda <- offset
  }
  if (!all(is.finite(log.lambda))) stop("log.lambda is not finite")
  log.lambda
}
