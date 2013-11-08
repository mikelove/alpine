# for imported BigWig files, this function adds GRanges with score 0
# before and after every contiguous region, so that plotting doesn't
# generate large diagonals
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
  gene <- sort(gene)
  z <- lapply(seq_along(gene), function(i) {
    shiftedGR <- shift(x[x %over% gene[i]], 
                       -start(gene[i]) + 
                         ifelse(i==1,0,sum(width(gene[seq_len(i-1)]))))
    if (addZero) {
      return(c(shiftedGR, GRanges(seqnames(x)[1],
                                  IRanges(sum(width(gene[seq_len(i)])),width=0),
                                  score=0)))
    } else {
      return(shiftedGR) 
    }
    }
  })
  do.call(c,z)
}
