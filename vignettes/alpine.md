<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{alpine}
-->

# Modeling and correcting fragment sequence bias

Here we show a brief example of using the *alpine* package to model
bias parameters and then using those parameters to estimate transcript
abundance. We load a subset of reads from four samples from the
GEUVADIS project. For more details on these files, see
`?ERR188297` in the *alpineData* package.




```r
library(alpine)
dir <- system.file("inst/extdata",package="alpine")
metadata <- read.csv(file.path(dir,"metadata.csv"),
                     stringsAsFactors=FALSE)
bam.files <- file.path(dir,paste0(metadata$Title,"_galignpairs.bam"))
names(bam.files) <- metadata$Title
stopifnot(all(file.exists(bam.files)))
```

To fit the bias model, we need to identify single-isoform genes.
We used the following chunk of code (here not
evaluated to generate a *GRangesList* of exons per single-isoform gene.


```r
library(ensembldb)
gtf.file <- "~/proj/alpine/alpine/inst/extdata/Homo_sapiens.GRCh38.84.gtf"
#gtf.file <- "Homo_sapiens.GRCh38.84.gtf"
txdb <- EnsDb(gtf.file) # already an EnsDb
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]
# pre-selected genes based on medium to high counts
# calculated using Rsubread::featureCounts
selected.genes <- scan("~/proj/alpine/alpine/inst/extdata/selected.genes.txt",what="char")
one.iso.txs <- txdf$tx_id[txdf$gene_id %in%
                          intersect(one.iso.genes, selected.genes)]
ebt0 <- exonsBy(txdb, by="tx")
ebt <- ebt0[one.iso.txs]
#save(ebt, file="ebt.rda")
```

Here we pick a subset of single-isoform genes based on the
number of exons, and the length. We show in comments the recommended
parameters to use in selecting this subset of genes,
although here we use different parameters. This choice of
different parameters here is to ensure the building of the
vignette takes a short period of time and does not use
much memory.


```r
library(GenomicRanges)
dir <- system.file("data",package="alpine")
load(file.path(dir,"ebt.rda"))
# more than 1 exon
table(elementNROWS(ebt))
```

```
## 
##  1  2  3  4  5  6  7  8  9 10 11 12 13 16 17 22 23 
## 61 55 17  8 15  9 11  4  4  5  2  3  2  1  1  1  1
```

```r
ebt <- ebt[elementNROWS(ebt) > 1]
# filter small genes and long genes
min.bp <- 1000 # better 800 bp
max.bp <- 2000 # better 5000 bp
gene.lengths <- sum(width(ebt))
summary(gene.lengths)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     362    1401    2967    3219    4715    9920
```

```r
ebt <- ebt[gene.lengths > min.bp & gene.lengths < max.bp]
length(ebt)
```

```
## [1] 22
```

```r
set.seed(1)
ebt <- ebt[sample(length(ebt),10)] # better 100 genes
```

# Fitting the bias model

Robust fitting of these bias
parameters is best with ~100 medium to highly expressed genes.
It is required to specify a minimum and maximum fragment
size (some lower and upper quantiles of the fragment length
distribution) and the read length. `minsize` and `maxsize`
are recommended to be the 2.5% and 97.5% of the fragment
length distribution. Currently
*alpine* only supports unstranded, paired-end RNA-seq with
fixed read length. Differences of +/- 1 bp in read length
across samples can be ignored.


```r
library(alpine)
library(BSgenome.Hsapiens.NCBI.GRCh38)
minsize <- 100 # better 80 for this data
maxsize <- 250 # better 350 for this data
readlength <- 75 
gene.names <- names(ebt)
names(gene.names) <- gene.names
```

The following function builds a list of information
about the fragment types from each gene in our
training set.


```r
system.time({
fragtypes <- lapply(gene.names, function(gene.name) {
                      buildFragtypes(exons=ebt[[gene.name]],
                                     genome=Hsapiens,
                                     readlength=readlength,
                                     minsize=minsize,
                                     maxsize=maxsize)
                    })
})
```

```
##    user  system elapsed 
##  16.506   2.485  19.275
```

```r
object.size(fragtypes)/1e6
```

```
## 243.367448 bytes
```

We can examine the information for a single gene:


```r
head(fragtypes[[1]], 3)
```

```
## DataFrame with 3 rows and 18 columns
##       start       end     relpos   fraglen        id fivep.test
##   <integer> <integer>  <numeric> <integer> <IRanges>  <logical>
## 1         1       100 0.02676660       100  [1, 100]      FALSE
## 2         1       101 0.02730193       101  [1, 101]      FALSE
## 3         1       102 0.02730193       102  [1, 102]      FALSE
##            fivep threep.test                threep        gc   GC40.90
##   <DNAStringSet>   <logical>        <DNAStringSet> <numeric> <numeric>
## 1  AGTTTCGGAACCC        TRUE GTGCCCGGGCTGCAGAAGCCC 0.6600000         0
## 2  AGTTTCGGAACCC        TRUE GGTGCCCGGGCTGCAGAAGCC 0.6633663         0
## 3  AGTTTCGGAACCC        TRUE GGGTGCCCGGGCTGCAGAAGC 0.6666667         0
##     GC40.80   GC20.90   GC20.80    gstart      gend gread1end gread2start
##   <numeric> <numeric> <numeric> <integer> <integer> <integer>   <integer>
## 1         0         0         1  33567125  33567026  33567051    33567100
## 2         0         0         1  33567125  33567025  33567051    33567099
## 3         0         0         1  33567125  33567024  33567051    33567098
```

# Defining bias models

The definition of bias models is extremely flexible in *alpine*.  The
`models` argument should be given as a list, where each element is
model.  The model itself should be provided as a list with elements
`formula` and `offset`. `offset` can be set to `NULL`.  The allowable
offsets are `fraglen` and/or `vlmm` which should be listed as a
character vector.

Any kind of R formula can be provided to formula, making use of the
fragment features, `gc` (0 to 1), `relpos` (0 to 1), `GC40.80`,
`GC40.90`, `GC20.80`, `GC20.90` (all indicator variables),
all stored in the elements contained in
`fragtypes`. Interactions between these terms and offsets,
e.g. `gc:fraglen` can also be fit.
We recommend providing formula as character vectors,
which are converted internally into formula, due to details in how R
formula make copies of objects from the environment.


```r
models <- list(
  "all" = list(formula = "count ~
  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  gene",
  offset=c("fraglen","vlmm"))
)
```

Here we fit a bias model
using fragment length, read start (VLMM),
fragment GC content, GC runs, relative position, and a term for
differences in expression across the genes (`+ gene`).
The knots and boundary knots for GC content
and relative position splines are fixed internally.
The returned object, `fitpar`, stores the information
as a list of fitted parameters across samples.


```r
system.time({
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
})
```

```
##    user  system elapsed 
##  71.947   6.075  79.162
```

```r
#save(fitpar, file="fitpar.rda")
```

# Visually exploring the bias parameters

Note that with more basepairs between `minsize` and `maxsize` and with
more genes used for estimation, the bias parameters 
would be more precise. As estimated here, they likely have high
variance from too few observations (paired-end fragments) across too
few genes.

First we set a palette to distinguish between samples


```r
library(RColorBrewer)
palette(brewer.pal(8,"Dark2"))
```

The fragment length distribution:


```r
perf <- as.integer(factor(metadata$Performer))
plotFragLen(fitpar, col=perf)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

The fragment GC bias curves:


```r
plotGC(fitpar, model="all", col=perf)
```

![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png)

The relative positional bias:


```r
plotRelPos(fitpar, model="all", col=perf)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

A 0-order version of the VLMM (note that the
VLMM that is used in the model includes positions
that are 1- and 2-order, so this plot does not
represent the final VLMM).


```r
plotOrder0(fitpar[["ERR188297"]][["vlmm.fivep"]][["order0"]])
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

```r
plotOrder0(fitpar[["ERR188297"]][["vlmm.threep"]][["order0"]])
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-2.png)


```r
print(head(fitpar[["ERR188297"]][["summary"]][["all"]]), row.names=FALSE)
```

```
##    Estimate Std. Error z value Pr(>|z|)
##  -3.8906360  1.0396561 -3.7422 1.82e-04
##   0.0306145  1.0036296  0.0305 9.76e-01
##  -0.4398308  0.6376371 -0.6898 4.90e-01
##  -3.0582911  2.0660590 -1.4803 1.39e-01
##  -8.2226083  0.8347416 -9.8505 6.82e-23
##   1.0157609  0.1626099  6.2466 4.19e-10
```

# Estimating transcript abundances

We pick a subset of genes for estimating transcript abundances.  If
the gene annotation includes genes with transcripts which span
multiple chromosomes or which do not have any overlap and are very far
apart, `splitGenesAcrossChroms` and `splitLongGenes`, respectively,
can be used to split these.  For again merging any overlapping
transcripts into "genes", the `mergeGenes` function can be used.  Here
we use the ENSEMBL gene annotation as is.


```r
one.iso.genes <- intersect(names(tab)[tab == 1], selected.genes)
two.iso.genes <- intersect(names(tab)[tab == 2], selected.genes)
three.iso.genes <- intersect(names(tab)[tab == 3], selected.genes)
set.seed(1)
subset.genes <- c(sample(one.iso.genes, 2),
                  sample(two.iso.genes, 2),
                  sample(three.iso.genes, 2)) 
txdf.sub <- txdf[txdf$gene_id %in% subset.genes,]
ebt.sub <- ebt0[txdf.sub$tx_id]
#save(ebt.sub, txdf.sub, subset.genes, file="estimationObjects.rda")
```


```r
dir <- system.file("data",package="alpine")
load(file.path(dir,"estimationObjects.rda"))
```

We specify a set of models. Any formula used here must
be equal to ones used in a fitted model of the same
name, except `+ gene` is replaced with `+ 0`.


```r
models <- list(
  "null"=list(formula=NULL, offset=NULL),
  "all"=list(formula="count~
  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  0",
  offset=c("fraglen","vlmm"))
  )
```

Here we estimate FPKM-scale abundances for multiple genes and multiple
samples.


```r
system.time({
res <- lapply(subset.genes, function(gene.name) {
         txs <- txdf.sub$tx_id[txdf.sub$gene_id == gene.name]
         estimateTheta(transcripts=ebt.sub[txs],
                       bamfiles=bam.files,
                       fitpar=fitpar,
                       genome=Hsapiens,
                       models=models,
                       readlength=readlength,
                       minsize=minsize,
                       maxsize=maxsize)
       })
})
```

```
##    user  system elapsed 
## 133.787   9.507 144.049
```

Each element of this list has the results for a single gene across all
samples, all models, and all isoforms of the gene:


```r
res[[1]][["ERR188297"]][["all"]]
```

```
## $theta
## ENST00000272233 
##        8031.435 
## 
## $lambda
## ENST00000272233 
##      0.02098951
```

```r
res[[6]][["ERR188297"]][["all"]]
```

```
## $theta
## ENST00000422252 ENST00000331479 ENST00000484235 
##       11346.488        8434.599       28801.224 
## 
## $lambda
## ENST00000422252 ENST00000331479 ENST00000484235 
##      0.03337403      0.02488118      0.02529331
```

These estimates use a default library size of 1e6, but can
be collated and rescaled using the total number of fragments
across all genes in `res` using the `extractAlpine` function.
This also will scale the estimates such that the total bias observed
over all genes is centered at 1 (log bias equal to zero).


```r
mat <- extractAlpine(res,
                     model="all",
                     nsamp=length(bam.files))
mat
```

```
##                    ERR188297    ERR188088    ERR188204    ERR188317
## ENST00000272233 4.408219e+04 9.764919e+04 5.823250e+04 5.485479e+04
## ENST00000374518 4.300741e+04 3.858410e+04 4.290873e+04 4.742574e+04
## ENST00000373758 3.123748e+04 2.805677e+04 2.894547e+04 3.306017e+04
## ENST00000634963 3.955291e-06 2.356548e+03 1.435929e-14 1.003936e-06
## ENST00000416247 8.893390e+03 6.891867e+03 1.028256e+04 9.450473e+03
## ENST00000376935 5.445548e+03 5.235020e+03 2.800799e+03 8.851717e+03
## ENST00000613913 7.015160e-04 3.271283e-04 4.402587e+02 4.326308e+02
## ENST00000353704 7.589523e+04 9.123782e+04 7.261348e+04 6.812086e+04
## ENST00000486056 8.040640e+03 5.844470e+03 1.882918e+01 2.552673e+03
## ENST00000422252 6.227754e+04 8.290126e+04 1.899833e+05 1.184271e+05
## ENST00000331479 4.629504e+04 2.709261e+04 2.945727e+04 2.707870e+04
## ENST00000484235 1.580815e+05 1.651025e+05 1.109300e+05 1.093027e+05
```

If we provide a *GRangesList* which contains the exons for each
transcript, the returned object will be a *SummarizedExperiment*.
The *GRangesList* does not have to be in the correct order,
the transcripts will be extracted by name to match the rows of the
FPKM matrix.


```r
se <- extractAlpine(res,
                    model="all",
                    nsamp=length(bam.files),
                    transcripts=ebt.sub)
se
```

```
## class: RangedSummarizedExperiment 
## dim: 12 4 
## metadata(0):
## assays(1): FPKM
## rownames(12): ENST00000272233 ENST00000374518 ... ENST00000331479
##   ENST00000484235
## rowData names(0):
## colnames(4): ERR188297 ERR188088 ERR188204 ERR188317
## colData names(0):
```

The matrix of FPKM values can be scaled using the median ratio method
of DESeq with the `normalizeDESeq` function. This is a robust method
which removes systematic differences in values across samples, and is
more appropriate than using the total count which is highly sensitive
to large abundance estimates from individual transcripts.


```r
norm.mat <- normalizeDESeq(mat, cutoff=0.1)
```


```r
sessionInfo()
```

```
## R Under development (unstable) (2016-03-21 r70361)
## Platform: x86_64-apple-darwin14.5.0 (64-bit)
## Running under: OS X 10.10.5 (Yosemite)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices datasets  utils    
## [8] methods   base     
## 
## other attached packages:
##  [1] alpine_0.1.5                          
##  [2] RColorBrewer_1.1-2                    
##  [3] GenomicAlignments_1.9.4               
##  [4] Rsamtools_1.25.0                      
##  [5] BSgenome.Hsapiens.NCBI.GRCh38_1.3.1000
##  [6] BSgenome_1.41.2                       
##  [7] rtracklayer_1.33.7                    
##  [8] Biostrings_2.41.4                     
##  [9] XVector_0.13.2                        
## [10] SummarizedExperiment_1.3.5            
## [11] ensembldb_1.5.8                       
## [12] GenomicFeatures_1.25.14               
## [13] AnnotationDbi_1.35.3                  
## [14] Biobase_2.33.0                        
## [15] GenomicRanges_1.25.8                  
## [16] GenomeInfoDb_1.9.1                    
## [17] IRanges_2.7.11                        
## [18] S4Vectors_0.11.5                      
## [19] BiocGenerics_0.19.1                   
## [20] magrittr_1.5                          
## [21] knitr_1.13                            
## [22] testthat_1.0.2                        
## [23] devtools_1.11.1                       
## [24] BiocInstaller_1.23.6                  
## 
## loaded via a namespace (and not attached):
##  [1] splines_3.4.0                 lattice_0.20-33              
##  [3] htmltools_0.3.5               interactiveDisplayBase_1.11.3
##  [5] XML_3.98-1.4                  RBGL_1.49.1                  
##  [7] withr_1.0.2                   DBI_0.4-1                    
##  [9] BiocParallel_1.7.4            speedglm_0.3-1               
## [11] stringr_1.0.0                 zlibbioc_1.19.0              
## [13] memoise_1.0.0                 evaluate_0.9                 
## [15] biomaRt_2.29.2                httpuv_1.3.3                 
## [17] markdown_0.7.7                Rcpp_0.12.5                  
## [19] xtable_1.8-2                  formatR_1.4                  
## [21] graph_1.51.0                  mime_0.4                     
## [23] AnnotationHub_2.5.4           digest_0.6.9                 
## [25] stringi_1.1.1                 shiny_0.13.2                 
## [27] grid_3.4.0                    tools_3.4.0                  
## [29] bitops_1.0-6                  RCurl_1.95-4.8               
## [31] RSQLite_1.0.0                 crayon_1.3.1                 
## [33] MASS_7.3-45                   Matrix_1.2-6                 
## [35] httr_1.2.0                    roxygen2_5.0.1               
## [37] R6_2.1.2                      compiler_3.4.0
```


