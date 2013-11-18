gr <- GRanges("chr1",IRanges(15e6+1,width=5e5))
fls <- list.files(".","*.bw")

getChunk <- function(gr, fls) {
  cvrList <- lapply(fls, function(f) {
    bwf <- BigWigFile(f)
    x <- import(bwf,which=gr)
    cvr <- coverage(x,weight=mcols(x)$score)[[as.character(seqnames(gr))]]
    as.integer(Views(cvr, ranges(gr))[[1]])
  })
  do.call(cbind, cvrList)
}

z <- getChunk(gr, fls)