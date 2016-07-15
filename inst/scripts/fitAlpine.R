metadata <- read.csv("../extdata/metadata.csv", stringsAsFactors=FALSE)
bam.files <- paste0("out/",metadata$Title,"_galignpairs.bam")
names(bam.files) <- metadata$Title
stopifnot(all(file.exists(bam.files)))

if (FALSE) {

library(ensembldb)
gtf.file <- "Homo_sapiens.GRCh38.84.gtf"
txdb <- EnsDb(basename(gtf.file))

txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]

# pre-selected genes
selected.genes <- scan("../extdata/selected.genes.txt",what="char")

one.iso.txs <- txdf$tx_id[txdf$gene_id %in% intersect(one.iso.genes, selected.genes)]
length(one.iso.txs)

ebt0 <- exonsBy(txdb, by="tx")
ebt <- ebt0[one.iso.txs]

save(ebt, file="ebt.rda")

}

library(GenomicRanges)
load("ebt.rda")

# more than 1 exon
table(elementNROWS(ebt))
ebt <- ebt[elementNROWS(ebt) > 1]

# filter small genes and long genes
min.bp <- 1000 # 800 bp
max.bp <- 2000 # 5000 bp
gene.lengths <- sum(width(ebt))
summary(gene.lengths)
ebt <- ebt[gene.lengths > min.bp & gene.lengths < max.bp]
length(ebt)

set.seed(1)
ebt <- ebt[sample(length(ebt),4)]

models <- list(
  "GC" = list(formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  gene",
  offset=c("fraglen"))
)

library(alpine)
library(BSgenome.Hsapiens.NCBI.GRCh38)

minsize <- 125 # better 80
maxsize <- 175 # better 350
readlength <- 75 

gene.names <- names(ebt)
names(gene.names) <- gene.names
fragtypes <- lapply(gene.names, function(gene.name) {
                     buildFragtypes(exons=ebt[[gene.name]],
                                    genome=Hsapiens,
                                    readlength=readlength,
                                    minsize=minsize,
                                    maxsize=maxsize,
                                    gc.str=FALSE,
                                    vlmm=FALSE)
                   })

fitpar <- lapply(bam.files, function(bf) {
                   fitBiasModels(genes=ebt,
                                 bamfile=bf,
                                 fragtypes=fragtypes,
                                 genome=Hsapiens,
                                 models=models,
                                 readlength=readlength,
                                 minsize=minsize,
                                 maxsize=maxsize)
                 })


save(fitpar, file="fitpar.rda")

sessionInfo()
