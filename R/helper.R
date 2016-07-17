#' Extract results from estimateTheta run across genes
#'
#' This function extracts estimates for a given model from a list
#' over many genes, returning a matrix with dimensions:
#' number of transcript x number of samples.
#' Here, the count of compatible fragments aligning to the
#' genes is used to estimate the FPKM, dividing out the previously
#' used estimate \code{lib.sizes}.
#' 
#' @param res a list where each element is the output of \link{estimateTheta}
#' @param model the name of a model, corresponds to names of \code{models}
#' used in \link{fitBiasModels}
#' @param nsamp the number of samples, corresponds to length of \code{bamfiles}
#' in \link{estimateTheta}
#' @param lib.sizes the vector of library sizes passed to \link{estimateTheta}.
#' not needed if \code{divideOut=FALSE}
#' @param divideOut logical, whether to divide out the initial estimate of
#' library size and to instead use the count of compatible fragments for
#' genes calculated by \link{estimateTheta}. Default is TRUE
#' @param transcripts an optional \code{GRangesList} of the exons for each
#' transcript. If this is provided, the output will be a
#' \code{SummarizedExperiment}. The transcripts do not need
#' to be provided in the correct order, \code{extractAlpine} will
#' find the correct transcript by the names in \code{res} and
#' put them in the correct order.
#' 
#' @return a matrix of FPKM values across transcripts and samples
#'
#' @export
extractAlpine <- function(res, model, nsamp, lib.sizes=1e6,
                          divideOut=TRUE, transcripts=NULL) {
  fpkm <- extractRes(res, model, "theta", nsamp)
  lambda <- extractRes(res, model, "lambda", nsamp)
  count <- extractRes(res, model, "count", nsamp)
  lambdaBar <- colMeans(lambda, na.rm=TRUE)
  colSumsCount <- colSums(count)
  multFactor <- if (divideOut) {
    lambdaBar * lib.sizes / colSumsCount
  } else {
    lambdaBar
  }
  mat <- sweep(fpkm, 2, multFactor, `*`)
  if (is.null(transcripts)) {
    return(mat)
  } else {
    row.ranges <- transcripts[rownames(mat)]
    se <- SummarizedExperiment(assays=list(FPKM=mat),
                               rowRanges=row.ranges)
    return(se)
  }
}

#' Split genes that have isoforms across chromosomes
#'
#' This function simply splits apart genes which have isoforms across multiple
#' chromosomes. New "genes" are created with the suffix "_cs" and a number.
#' 
#' @param ebg an exons-by-genes GRangesList, created with \code{exonsBy}
#' @param txdf a data.frame created by running \code{select} on a TxDb object.
#' Must have columns TXCHROM and GENEID
#'
#' @return a list of manipulated \code{ebg} and \code{txdf}
#'
#' @export
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

#' Split very long genes
#'
#' This function splits genes which have a very long range (e.g. 1 Mb),
#' and new "genes" are formed where each isoform is its own "gene",
#' with the suffix "_ls" and a number.
#' It makes sense to turn each isoform into its own gene only if this
#' function is followed by \link{mergeGenes}.
#' 
#' @param ebg an exons-by-genes GRangesList, created with \code{exonsBy}
#' @param ebt an exons-by-tx GRangesList, created with \code{exonsBy}
#' @param txdf a data.frame created by running \code{select} on a TxDb object.
#' Must have columns GENEID and TXNAME, where TXNAME corresponds to the
#' names of \code{ebt}. Note: this requires renaming \code{ebt}.
#' @param long a numeric value such that ranges longer than this are "long"
#'
#' @return a list of manipulated \code{ebg} and \code{txdf}
#'
#' @export
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

#' Merge overlapping "genes" into gene clusters
#'
#' This function looks for overlapping exons in \code{ebg}.
#' The overlapping "genes" are used to form a graph.
#' Any connected components in the graph (sets of "genes"
#' which can be reached from each other through overlap relations)
#' are connected into a new gene cluster, which is given the
#' suffix "_mrg" and using one of the original gene names.
#' 
#' @param ebg an exons-by-genes GRangesList, created with \code{exonsBy}
#' @param txdf a data.frame created by running \code{select} on a TxDb object.
#' Must have a column GENEID.
#' @param ignore.strand Default is TRUE.
#'
#' @return a manipulated \code{txdf}.
#'
#' @export
mergeGenes <- function(ebg, txdf, ignore.strand=TRUE) {
  fo <- findOverlaps(ebg, ignore.strand=ignore.strand)
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

#' DESeq median ratio normalization for matrix
#'
#' Simple implementation of DESeq median ratio normalization
#'
#' @param mat a matrix of numeric values
#' @param cutoff a numeric value to be used as the cutoff
#' for the row means of \code{mat}. Only rows with row mean
#' larger than \code{cutoff} are used for calculating
#' the size factors
#'
#' @return a matrix with the median ratio size factors
#' divided out
#'
#' @references Anders, S. and Huber, W.,
#' Differential expression analysis for sequence count data.
#' Genome Biology (2010) doi: 10.1186/gb-2010-11-10-r106
#' 
#' @export
normalizeDESeq <- function(mat, cutoff) {
  mat2 <- mat[rowMeans(mat) > cutoff,,drop=FALSE]
  loggeomeans <- rowMeans(log(mat2))
  logratio <- (log(mat2) - loggeomeans)[is.finite(loggeomeans),,drop=FALSE]
  sf <- exp(apply(logratio, 2, median, na.rm=TRUE))
  sweep(mat, 2, sf, "/")
}

#' Get fragment widths
#'
#' From a BAM file and a particular transcript (recommened
#' to be the single isoform of a gene), this function
#' returns estimates of the fragment widths, by mapping the
#' fragment alignments to the transcript coordinates.
#'
#' @param bamfile a character string pointing to a BAM file
#' @param tx a GRanges object of the exons of a single isoform gene
#'
#' @return a numeric vector of estimated fragment widths
#'
#' @export
getFragmentWidths <- function(bamfile, tx) {
  gap <- readGAlignmentPairs(bamfile, param=ScanBamParam(which=range(tx)))
  stopifnot(length(gap) > 0)
  fo <- findCompatibleOverlaps(gap, GRangesList(tx=tx))
  stopifnot(length(fo) > 0)
  gap <- gap[queryHits(fo)]
  left <- first(gap)
  right <- last(gap)
  first.minus <- as.vector(strand(first(gap)) == "-")
  left[first.minus] <- last(gap)[first.minus]
  right[first.minus] <- first(gap)[first.minus]
  left.tx <- start(mapToTranscripts(GRanges(seqnames(gap),
                                            IRanges(start(left),width=1)),
                                    GRangesList(tx=tx)))
  right.tx <- end(mapToTranscripts(GRanges(seqnames(gap),
                                           IRanges(end(right),width=1)),
                                   GRangesList(tx=tx)))
  w <- right.tx - left.tx + 1
  if (as.character(strand(tx)[1]) == "-") {
    w <- w * -1
  }
  return(w)
}

######### unexported helper functions #########

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
getReadLength <- function(bamfiles) {
  getRL1 <- function(file) {
    qwidth(readGAlignments(BamFile(file, yieldSize=1)))
  }
  sapply(bamfiles, getRL1)
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
