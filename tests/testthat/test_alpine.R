context("alpine")
test_that("alpine works", {
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

  
  # model missing '+ gene' gives error
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

  
  # works to fit only fraglen and vlmm
  models <- list("readstart"=list(formula=NULL,offset=c("fraglen","vlmm")))
  fitpar <- fitBiasModels(
    genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
    models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
  )

  
  # works to fit different knots
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

  
  # works to estimate abundances
  models <- list(
    "GC"=list(formula="count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + gene",
              offset=c("fraglen","vlmm"))
  )
  fitpar <- list(
    ERR188088 = fitBiasModels(
      genes=ebt.fit[gene.names],bam.file=bam.file,fragtypes=fragtypes,genome=Hsapiens,
      models=models,readlength=readlength,minsize=minsize,maxsize=maxsize
    )
  )
  # specify models using a character vector 
  model.names <- c("fraglen","fraglen.vlmm","GC")
  txs <- txdf.theta$tx_id[txdf.theta$gene_id == "ENSG00000198918"]
  res <- estimateAbundance(transcripts=ebt.theta[txs],
                           bam.files=bam.file,
                           fitpar=fitpar,
                           genome=Hsapiens,
                           model.names=model.names)

  
  # works to predict coverage
  pred.cov <- predictCoverage(gene=ebt.fit[["ENST00000379660"]],
                              bam.files=bam.file,
                              fitpar=fitpar.small,
                              genome=Hsapiens,
                              model.names=model.names)
  # plot
  frag.cov <- pred.cov[["ERR188088"]][["frag.cov"]]
  plot(frag.cov, type="l", lwd=3, ylim=c(0,max(frag.cov)*1.5))
  for (i in seq_along(model.names)) {
    m <- model.names[i]
    pred <- pred.cov[["ERR188088"]][["pred.cov"]][[m]]
    lines(pred/mean(pred)*mean(frag.cov), col=i+1, lwd=3)
  }
  legend("topright", legend=c("observed",model.names),
         col=seq_len(length(model.names)+1), lwd=3)
})
