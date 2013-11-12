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

# this function takes a GRanges x with coverage info and a
# GRanges gene with exons, and moves the coverage to start
# at the origin, and removes introns
exonCoverage <- function(x, gene, addZero=TRUE) {
  x <- sort(x)
  gene <- reduce(sort(gene))
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

# takes a GRanges x and a GRanges gene
# returns a data.frame of base by base coverage along exons
# from TSS continuing downstream, i.e. negative strand genes 
# have the wiggle flipped
GRangesCoverageToNumericCoverage <- function(x, gene) {
  if (length(x) == 0) return(data.frame(x=numeric(0),y=numeric(0)))
  totalWidth <- sum(width(gene))
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
