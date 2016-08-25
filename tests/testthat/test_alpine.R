library(alpineData)
library(GenomicAlignments)
library(rtracklayer)
gap <- ERR188088()
dir <- system.file(package="alpineData", "extdata")
bam.file <- c("ERR188088" = file.path(dir,"ERR188088.bam"))
export(gap, con=bam.file)
library(GenomicRanges)
library(BSgenome.Hsapiens.NCBI.GRCh38)
data(preprocessedData)
readlength <- 75
minsize <- 125
maxsize <- 175
gene.names <- names(ebt.fit)[6:8]
names(gene.names) <- gene.names
fragtypes <- lapply(gene.names, function(gene.name) {
                      buildFragtypes(ebt.fit[[gene.name]],
                                     Hsapiens, readlength,
                                     minsize, maxsize)
})

# check model missing '+ gene'
models <- list(
  "GC"=list(formula="count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk)",
            offset=c("fraglen","vlmm"))
)
expect_error(
  fitBiasModels(
    genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
    models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
  )
)

# check fitting only fraglen and vlmm
models <- list("readstart"=list(formula=NULL,offset=c("fraglen","vlmm")))
fitpar <- fitBiasModels(
  genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
  models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
)

# check fitting different knots
models <- list(
  "GC"=list(formula="count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + gene",
            offset=c("fraglen","vlmm"))
)
fitpar <- fitBiasModels(
  genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
  models=models,readlength=readlength,minsize=minsize,maxsize=maxsize,
  gc.knots=seq(from=.3,to=.6,length=5), gc.bk=c(0,1)
)
plotGC(fitpar, model="GC")

# test estimate abundances
models <- list(
  "GC"=list(formula="count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + gene",
            offset=c("fraglen","vlmm"))
)
fitpar <- fitBiasModels(
  genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
  models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
)
models <- list("fraglen"="fraglen",
               "readstart"=c("fraglen","vlmm"),
               "GC"="GC")
## res <- estimateAbundance(transcripts=ebt.theta[txs],
##                          bam.files=bam.file,
##                          fitpar=fitpar.small,
##                          genome=Hsapiens,
##                          models=models)
