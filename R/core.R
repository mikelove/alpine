# for imported BigWig files, this function adds GRanges with score 0
# before and after every contiguous region, so that plotting doesn't
# generate diagonal lines over large regions
addGRangesAtZero <- function(x) {
  x <- sort(x)
  idx <- which(start(x)[-1] != end(x)[-length(x)] + 1)
  zeros <- GRanges(seqnames(x)[1],
                   IRanges(c(end(x)[idx] + 1,start(x)[idx+1] - 1),width=0),
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
                    IRanges(cumsum(width(gene))[-length(gene)],width=0),
                    score=rep(0,length(gene)-1))
    return(sort(c(z,zero)))
  } 
}

# takes GRanges and returns data.frame of base by base coverage
GRangesCoverageToNumericCoverage <- function(x, gene) {
  if (length(x) == 0) return(data.frame(x=numeric(0),y=numeric(0)))
  totalWidth <- sum(width(gene))
  strand <- as.character(strand(gene)[1])
  d0 <- data.frame(x=start(x), y=mcols(x)$score)
  lastPos <- d0$x[length(d0$x)]
  widths <- c(d0$x, lastPos+1, totalWidth+1) - c(0, d0$x, lastPos)
  rlex <- Rle(values=c(0,d0$y,0), lengths=widths)
  if (strand == "-") {
    y <- rev(as.numeric(rlex))
  } else {
    y <- as.numeric(rlex)
  }
  data.frame(x=seq_len(totalWidth-1), y=y)
}
