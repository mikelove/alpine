# for all sequences (e.g. chromsomes)
# bwGR is a GRanges produced by import(BigWig(...))
# exons is a GRangesList of exons by genes produced by exonsBy(txdb, by="gene")
# this should also be sorted using exonsSort <- endoapply(exons, sort)
# this returns a list of lists:
# for each sequence for each gene an integer vector of exonic coverage
bwGRInExonsGRListToInt <- function(bwGR, exons, flipNames=NULL) {
  x <- coverage(bwGR,weight=mcols(bwGR)$score)
  seqNames <- names(x)
  names(seqNames) <- seqNames
  lapply(seqNames, function(seqname) {
    exonsSub <- exons[seqnames(exons) == seqname]
    exonsSub <- exonsSub[sapply(exonsSub,length) > 0]
    intCov <- rleInRangesListToInt(x[[seqname]],ranges(exonsSub))
    if (!is.null(flipNames)) {
      idx <- names(intCov) %in% flipNames
      intCov[idx] <- lapply(intCov[idx], rev) 
    }
    intCov
  })
}

# for a single sequence (e.g. chromosome)
# x is an Rle
# y is a RangesList
rleInRangesListToInt <- function(x,y) {
  sp <- rep(names(y),sapply(y,length))
  ry <- unlist(y,use.names=FALSE)
  intCovList <- split(viewApply(Views(x,ry),as.integer),sp)
  lapply(intCovList, function(z) do.call(c,z))
}

# for imported BigWig files, this function adds GRanges with score 0
# before and after every contiguous region, so that plotting doesn't
# generate diagonal lines over large regions
addGRangesAtZero <- function(x) {
  x <- sort(x)
  idx <- which(start(x)[-1] != end(x)[-length(x)] + 1)
  zeros <- GRanges(seqnames(x)[1],
                   IRanges(c(end(x)[idx] + 1,start(x)[idx+1] - 1),width=1),
                 score=rep(0,2*length(idx)))
  sort(c(x,zeros))
}

# NOTE: slow, designed for plotting
# this function takes a GRanges x with coverage info and a
# GRanges gene with exons, and moves the coverage to start
# at the origin, and removes introns
exonCoverage <- function(x, gene, addZero=TRUE) {
  x <- sort(x)
  gene <- reduce(gene)
  fo <- findOverlaps(x,gene)
  xR <- restrict(x[queryHits(fo)], 
                 start=start(gene)[subjectHits(fo)],
                 end=end(gene)[subjectHits(fo)])
  foR <- findOverlaps(xR, gene)
  shiftAmount <- -start(gene) + c(0,cumsum(width(gene))[-length(gene)]) + 1
  z <- shift(xR[queryHits(foR)], shiftAmount[subjectHits(foR)])
  if (!addZero) {
    return(z)
  } else {
    zero <- GRanges(seqnames(x)[1],
                    IRanges(cumsum(width(gene))[-length(gene)],width=1),
                    score=rep(0,length(gene)-1))
    return(sort(c(z,zero)))
  } 
}

# NOTE: slow, designed for plotting
# takes a GRanges produced by exonCoverage x and a GRanges gene
# returns a numeric of base by base coverage along exons
# from TSS continuing downstream, i.e. negative strand genes 
# have the wiggle flipped
exonCoverageToNumericCoverage <- function(x, gene) {
  x <- sort(x)
  gene <- reduce(gene)
  totalWidth <- sum(width(gene))
  if (length(x) == 0) return(numeric(totalWidth))
  x <- restrict(x,start=1,end=totalWidth)
  strand <- as.character(strand(gene)[1])
  firstPos <- start(x)[1]
  lastPos <- end(x)[length(x)]
  if (firstPos == 1 & lastPos == totalWidth) {
    rlex <- Rle(values=mcols(x)$score, lengths=width(x))
  } else if (firstPos == 1) {
    widths <- c(width(x), totalWidth - lastPos)
    rlex <- Rle(values=c(mcols(x)$score, 0), lengths=widths)
  } else if (lastPos == totalWidth) {
    widths <- c(firstPos - 1, width(x))
    rlex <- Rle(values=c(0, mcols(x)$score), lengths=widths)
  } else {
    widths <- c(firstPos - 1, width(x), totalWidth - lastPos)
    rlex <- Rle(values=c(0, mcols(x)$score, 0), lengths=widths)
  }
  if (strand == "-") {
    return(rev(as.numeric(rlex)))
  } else {
    return(as.numeric(rlex))
  }
}


