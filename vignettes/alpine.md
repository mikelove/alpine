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
dir <- system.file("extdata",package="alpine")
metadata <- read.csv(file.path(dir,"metadata.csv"),
                     stringsAsFactors=FALSE)
bam.files <- file.path(dir,paste0(metadata$Title,"_galignpairs.bam"))
names(bam.files) <- metadata$Title
stopifnot(all(file.exists(bam.files)))
```

To fit the bias model, we need to identify single-isoform genes.
We used the following chunk of code (here not
evaluated) to generate a *GRangesList* of exons per single-isoform gene.


```r
library(ensembldb)
gtf.file <- "~/proj/alpine/Homo_sapiens.GRCh38.84.gtf"
#gtf.file <- "Homo_sapiens.GRCh38.84.gtf"
txdb <- EnsDb(gtf.file) # already an EnsDb
txdf <- transcripts(txdb, return.type="DataFrame")
tab <- table(txdf$gene_id)
one.iso.genes <- names(tab)[tab == 1]
# pre-selected genes based on medium to high counts
# calculated using Rsubread::featureCounts
selected.genes <- scan("~/proj/alpine/alpine/inst/extdata/selected.genes.txt",
                       what="char")
one.iso.txs <- txdf$tx_id[txdf$gene_id %in%
                          intersect(one.iso.genes, selected.genes)]
ebt0 <- exonsBy(txdb, by="tx")
ebt.fit <- ebt0[one.iso.txs]
#save(ebt.fit, file="~/proj/alpine/alpine/data/ebtfit.rda")
```

Here we pick a subset of single-isoform genes based on the
number of exons, and the length. We show in comments the recommended
parameters to use in selecting this subset of genes,
although here we use different parameters to ensure the building of
the vignette takes only a short period of time and does not use much memory.


```r
library(GenomicRanges)
```


```r
dir <- system.file("data",package="alpine")
load(file.path(dir,"ebtfit.rda"))
# filter small genes and long genes
min.bp <- 600
max.bp <- 7000 
gene.lengths <- sum(width(ebt.fit))
summary(gene.lengths)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    66.0   689.8  1585.0  2147.0  3352.0  6103.0
```

```r
ebt.fit <- ebt.fit[gene.lengths > min.bp & gene.lengths < max.bp]
length(ebt.fit)
```

```
## [1] 25
```

```r
set.seed(1)
# better to use ~100 genes
ebt.fit <- ebt.fit[sample(length(ebt.fit),10)] 
```

## Defining a set of fragment types

Robust fitting of these bias parameters is best with ~100 medium to
high count genes, e.g. mean count across samples between 200 and
10,000. These counts can be identified by *featureCounts* from the
*Rsubread* Bioconductor package, for example.
It is required to specify a minimum and maximum fragment size
which should be lower and upper quantiles of the fragment length
distribution. The `minsize` and `maxsize`
arguments are recommended to be roughly the 2.5% and 97.5% of the
fragment length distribution. This can be quickly estimated using the
helper function `getFragmentWidths`, iterating over a few
single-isoform genes with sufficient counts:


```r
w <- getFragmentWidths(bam.files[1], ebt.fit[[1]])
c(summary(w), Number=length(w))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  Number 
##    79.0   161.0   189.0   197.1   224.0   468.0  1762.0
```

```r
quantile(w, c(.025, .975))
```

```
##   2.5%  97.5% 
## 122.05 321.95
```

It is also required to specify the read length. Currently *alpine*
only supports unstranded, paired-end RNA-seq with fixed read
length. Differences of +/- 1 bp in read length across samples can be
ignored.


```r
getReadLength(bam.files)
```

```
## ERR188297 ERR188088 ERR188204 ERR188317 
##        75        75        76        76
```

Here we use a very limited range of fragment lengths for speed, but
for a real analysis we would suggest using the minimum and maximum
of the quantiles computed above across all samples (the minimum of the
lower quantiles and the maximum of the upper quantiles).


```r
library(alpine)
library(BSgenome.Hsapiens.NCBI.GRCh38)
```

```
## Loading required package: BSgenome
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## Loading required package: rtracklayer
```

```r
minsize <- 125 # better 80 for this data
maxsize <- 175 # better 350 for this data
readlength <- 75 
gene.names <- names(ebt.fit)
names(gene.names) <- gene.names
```

The following function builds a list of *DataFrames* which store
information about the fragment types from each gene in our
training set.


```r
system.time({
fragtypes <- lapply(gene.names, function(gene.name) {
                      buildFragtypes(exons=ebt.fit[[gene.name]],
                                     genome=Hsapiens,
                                     readlength=readlength,
                                     minsize=minsize,
                                     maxsize=maxsize,
                                     gc.str=FALSE)
                    })
})
```

```
##    user  system elapsed 
##  15.244   0.416  15.666
```

```r
object.size(fragtypes)/1e6
```

```
## 104.540624 bytes
```

We can examine the information for a single gene:


```r
head(fragtypes[[1]], 3)
```

```
## DataFrame with 3 rows and 14 columns
##       start       end     relpos   fraglen        id fivep.test
##   <integer> <integer>  <numeric> <integer> <IRanges>  <logical>
## 1         1       125 0.01032279       125  [1, 125]      FALSE
## 2         1       126 0.01032279       126  [1, 126]      FALSE
## 3         1       127 0.01048665       127  [1, 127]      FALSE
##            fivep threep.test                threep        gc    gstart
##   <DNAStringSet>   <logical>        <DNAStringSet> <numeric> <integer>
## 1  ATCCGGGGCAGCG        TRUE TCAATTCAAAATGTCTCAGGA 0.7680000  32277664
## 2  ATCCGGGGCAGCG        TRUE GTCAATTCAAAATGTCTCAGG 0.7619048  32277664
## 3  ATCCGGGGCAGCG        TRUE TGTCAATTCAAAATGTCTCAG 0.7637795  32277664
##        gend gread1end gread2start
##   <integer> <integer>   <integer>
## 1  32309735  32277738    32277714
## 2  32309736  32277738    32277715
## 3  32309737  32277738    32277716
```

## Defining and fitting bias models

The definition of bias models is extremely flexible in *alpine*.  The
`models` argument should be given as a list, where each element is
model.  The model itself should be provided as a list with elements
`formula` and `offset`. Either `formula` or `offset` can be set to
`NULL` for a given model. 
The allowable offsets are `fraglen` and/or `vlmm` which should be
provided in a character vector.
Offsets are only estimated once for all models, so setting
`formula=NULL` only makes sense if extra offsets are desired
which were not already calculated by other models.

Any kind of R formula can be provided to `formula`, making use of the
fragment features:

* `gc` (fragment GC content from 0 to 1)
* `relpos` (fragment midpoint relative position from 0 to 1)
* `GC40.80`, `GC40.90`, `GC20.80`, `GC20.90` (indicator variables
  indicating the presence of, e.g. a 40 bp stretch of 80% or higher GC
  content within the fragment)

These fragment features reference columns of information stored in
`fragtypes`.  Interactions between these terms and offsets are also
possible, e.g. `gc:fraglen`.  We recommend providing formula as
character vectors, which are converted internally into formula, due to
details in how R formula make copies of objects from the environment.


```r
models <- list(
  "GC" = list(formula = "count ~
  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  gene",
  offset=c("fraglen")),
  "all" = list(formula = "count ~
  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  gene",
  offset=c("fraglen","vlmm"))
)
```

Here we fit one bias model, `GC`, using fragment length, fragment GC
content, relative position,
and a term for differences in expression across the genes (`+
gene`).  We fit another bias model, `all`, with all the terms of the
first but additionally with read start bias (encoded by a Variable
Length Markov Model, or VLMM).  The knots and boundary knots for GC
content (`gc`) and relative position (`relpos`) splines are fixed
internally.  The returned object, `fitpar`, stores the information as
a list of fitted parameters across samples.


```r
system.time({
fitpar <- lapply(bam.files, function(bf) {
                   fitBiasModels(genes=ebt.fit,
                                 bam.file=bf,
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
##  53.860   0.792  54.784
```

```r
#fitpar.small <- fitpar # save as fitpar.small for examples
#save(fitpar.small, file="~/proj/alpine/alpine/data/fitparsmall.rda")
```

## Visually exploring the bias parameters

Note that with more basepairs between `minsize` and `maxsize` and with
more genes used for estimation, the bias parameters would be more
precise. As estimated here, they have high variance from too few
observations (paired-end fragments) across too few genes.

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

![plot of chunk fraglen](figure/fraglen-1.png)

The fragment GC bias curves:


```r
plotGC(fitpar, model="all", col=perf)
```

![plot of chunk gccurve](figure/gccurve-1.png)


The relative position curves:


```r
plotRelPos(fitpar, model="all", col=perf)
```

![plot of chunk relpos](figure/relpos-1.png)

A 0-order version of the VLMM (note that the VLMM that is used in the
model includes positions that are 1- and 2-order, so this plot does
not represent the final VLMM used in bias estimation or in estimation
of abundances).


```r
plotOrder0(fitpar[["ERR188297"]][["vlmm.fivep"]][["order0"]])
```

![plot of chunk vlmm](figure/vlmm-1.png)

```r
plotOrder0(fitpar[["ERR188297"]][["vlmm.threep"]][["order0"]])
```

![plot of chunk vlmm](figure/vlmm-2.png)

A coefficient table for the terms in `formula`:


```r
print(head(fitpar[["ERR188297"]][["summary"]][["all"]]), row.names=FALSE)
```

```
##    Estimate Std. Error z value Pr(>|z|)
##  -4.9780792  1.2944129 -3.8458 1.20e-04
##   1.3428276  1.2486577  1.0754 2.82e-01
##   1.8331308  0.8646923  2.1200 3.40e-02
##   1.8867834  2.5408684  0.7426 4.58e-01
##  -3.8256451  2.0065307 -1.9066 5.66e-02
##   0.8080533  0.1568161  5.1529 2.57e-07
```

## Estimating transcript abundances

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
genes.theta <- c(sample(one.iso.genes, 2),
                 sample(two.iso.genes, 2),
                 sample(three.iso.genes, 2)) 
txdf.theta <- txdf[txdf$gene_id %in% genes.theta,]
ebt.theta <- ebt0[txdf.theta$tx_id]
#save(ebt.theta, txdf.theta, genes.theta, file="~/proj/alpine/alpine/data/thetaobjs.rda")
```


```r
dir <- system.file("data",package="alpine")
load(file.path(dir,"thetaobjs.rda"))
```

We specify a set of models. Any formula used here must be equal to
those used in a fitted model of the same name, except `+ gene` is
replaced with `+ 0`.


```r
models <- list(
  "null"=list(formula=NULL, offset=NULL),
  "GC"=list(formula="count~
  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  0",
  offset=c("fraglen"))
  )
```

Here we estimate FPKM-scale abundances for multiple genes and multiple
samples. If `lib.sizes` is not specified, a default value of 1e6
is used. `estimateTheta` works one gene at a time, where the
`transcripts` argument expects a *GRangesList* of the exons for each
transcript (multiple if the gene has multiple isoforms).


```r
system.time({
res <- lapply(genes.theta, function(gene.name) {
         txs <- txdf.theta$tx_id[txdf.theta$gene_id == gene.name]
         estimateTheta(transcripts=ebt.theta[txs],
                       bam.files=bam.files,
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
##  44.836   0.316  45.226
```

```r
#save(res, file="~/proj/alpine/alpine/data/res.rda")
```

Each element of this list has the abundances (`theta`) and average
bias (`lambda`) for a single gene across all samples, all models, and all
isoforms of the gene: 


```r
res[[1]][["ERR188297"]][["GC"]]
```

```
## $theta
## ENST00000259030 
##       0.3425053 
## 
## $lambda
## ENST00000259030 
##        73.34301
```

```r
res[[6]][["ERR188297"]][["GC"]]
```

```
## $theta
## ENST00000477403 ENST00000468844 ENST00000361575 
##     7.869292934     0.008579681     3.014203873 
## 
## $lambda
## ENST00000477403 ENST00000468844 ENST00000361575 
##        82.26703        70.47448        82.04643
```

The `extractAlpine` function can be used to collate estimates from
across all genes.  `extractAlpine` will scale the estimates such that
the total bias observed over all transcripts is centered at 1.  The
estimates produce by `estimateTheta` presume a default library size of
1e6, but will be rescaled using the total number of fragments across
genes when using `extractAlpine` (if this library size rescaling is
not desired, choose `divide.out=FALSE`).


```r
mat <- extractAlpine(res, model="GC")
mat
```

```
##                    ERR188297   ERR188088   ERR188204    ERR188317
## ENST00000259030 5.557029e+03  17769.0617  10514.5581 1.302500e+04
## ENST00000304788 4.300769e+03   5620.1425   8585.7954 7.675728e+03
## ENST00000295025 1.341303e+04   9652.3313  18721.2112 1.531396e+04
## ENST00000394479 4.793712e+03   1000.5219   6205.6560 4.114905e+03
## ENST00000330871 1.525750e+04  33787.1234  10879.3160 1.139717e+04
## ENST00000587578 0.000000e+00      0.0000      0.0000 0.000000e+00
## ENST00000264254 3.832130e+04  50374.8018  53612.5766 6.396532e+04
## ENST00000416255 2.095944e+03   6057.3662   4859.2903 2.357083e+03
## ENST00000450127 5.805294e-32   1681.2494    585.0105 2.387932e-28
## ENST00000477403 1.276765e+05 128681.9577 243208.0857 1.707819e+05
## ENST00000468844 1.392023e+02    634.9873   1192.6475 1.167633e+03
## ENST00000361575 4.890440e+04  63053.9109  84102.2731 1.851644e+05
```

If we provide a *GRangesList* which contains the exons for each
transcript, the returned object will be a *SummarizedExperiment*.
The *GRangesList* provided to `transcripts` does not have to be in the
correct order, the transcripts will be extracted by name to match the
rows of the FPKM matrix.


```r
se <- extractAlpine(res, model="GC", transcripts=ebt.theta)
se
```

```
## class: RangedSummarizedExperiment 
## dim: 12 4 
## metadata(0):
## assays(1): FPKM
## rownames(12): ENST00000259030 ENST00000304788 ... ENST00000468844
##   ENST00000361575
## rowData names(0):
## colnames(4): ERR188297 ERR188088 ERR188204 ERR188317
## colData names(0):
```

The matrix of FPKM values can be scaled using the median ratio method
of DESeq with the `normalizeDESeq` function. This is a robust method
which removes systematic differences in values across samples, and is
more appropriate than using the total count which is sensitive to
very large abundance estimates for a minority of transcripts. 


```r
norm.mat <- normalizeDESeq(mat, cutoff=0.1)
```

## Simulating using empirically estimated GC bias

The fragment GC bias which *alpine* estimates can be used in
downstream simulations, for example in the *polyester* Bioconductor
package. All we need to do is to run the *plotGC* function, but
specifying that instead of a plot, we want to return a matrix of
probabilities for each percentile of fragment GC content. This matrix
can be provided to the `frag_GC_bias` argument of *simulate_experiment*.

We load a `fitpar` object that was run with the fragment length range
[80,350] bp. 


```r
data(fitpar)
prob.mat <- plotGC(fitpar, "all", return.type=2)
head(prob.mat)
```

```
##       ERR188297 ERR188088  ERR188204  ERR188317
## 0    0.04366855 0.2561645 0.06914584 0.07234787
## 0.01 0.04936226 0.2725952 0.07667866 0.08005226
## 0.02 0.05578805 0.2900464 0.08501986 0.08856473
## 0.03 0.06302707 0.3085437 0.09424127 0.09795502
## 0.04 0.07116603 0.3281071 0.10441771 0.10829555
## 0.05 0.08029675 0.3487502 0.11562636 0.11966079
```

If `return.type=0` (the default) the function makes the plot of log
fragment rate over fragment GC content. If `return.type=1` the
function returns the matrix of log fragment rate over percentiles of
fragment GC content, and if `return.type=2`, the matrix returns
probabilities of observing fragments based on percentiles of fragment
GC content (the log fragment rate exponentiated and scaled to have a
maximum of 1). The matrix returned by `return.type=2` is appropriate
for downstream use with *polyester*.

## Plotting predicted fragment coverage

In the *alpine* paper, it was shown that models incorporating fragment
GC bias can be a better predictor of test set RNA-seq fragment
coverage, compared to models incorporating read start bias. Here we
show how to predict fragment coverage for a single-isoform gene using
a variety of fitted bias models. The models involving formula need to
have the exact same name and form as a fitted model in `fitpar`.


```r
pred.models <- list(
  "fraglen" = list(formula=NULL, offset=c("fraglen")),
  "readstart" = list(formula=NULL, offset=c("fraglen","vlmm")),
  "GC" = list(formula = "count ~
  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  0",
  offset=c("fraglen")),
  "all" = list(formula = "count ~
  ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  0",
  offset=c("fraglen","vlmm"))
)
```

The following function computes the predicted coverage for one
single-isoform gene. We load a `fitpar` object that was run
with the fragment length range [80,350] bp.


```r
data(fitpar)
system.time({
  pred.cov <- predictCoverage(gene=ebt.fit[["ENST00000245479"]],
                              bam.files=bam.files["ERR188204"],
                              fitpar=fitpar,
                              genome=Hsapiens,
                              models=pred.models,
                              readlength=readlength,
                              minsize=80,
                              maxsize=350)
})
```

```
##    user  system elapsed 
##  22.828   0.288  23.206
```

We can plot the observed and predicted coverage for one of the
samples: 


```r
palette(brewer.pal(9, "Set1"))
frag.cov <- pred.cov[["ERR188204"]][["frag.cov"]]
plot(frag.cov, type="l", lwd=3, ylim=c(0,max(frag.cov)*1.5))
for (i in seq_along(pred.models)) {
  m <- names(pred.models)[i]
  pred <- pred.cov[["ERR188204"]][["pred.cov"]][[m]]
  lines(pred, col=i, lwd=3)
}
legend("topright", legend=c("observed",names(pred.models)),
       col=c("black",seq_along(pred.models)), lwd=3)
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)

## Session information


```r
sessionInfo()
```

```
## R Under development (unstable) (2016-05-23 r70660)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices datasets  utils    
## [8] methods   base     
## 
## other attached packages:
##  [1] RColorBrewer_1.1-2                    
##  [2] BSgenome.Hsapiens.NCBI.GRCh38_1.3.1000
##  [3] BSgenome_1.41.2                       
##  [4] rtracklayer_1.33.2                    
##  [5] Biostrings_2.41.1                     
##  [6] XVector_0.13.0                        
##  [7] GenomicRanges_1.25.0                  
##  [8] GenomeInfoDb_1.9.1                    
##  [9] IRanges_2.7.0                         
## [10] S4Vectors_0.11.1                      
## [11] BiocGenerics_0.19.0                   
## [12] alpine_0.1.6                          
## [13] magrittr_1.5                          
## [14] knitr_1.13                            
## [15] devtools_1.11.1                       
## [16] BiocInstaller_1.23.6                  
## 
## loaded via a namespace (and not attached):
##  [1] formatR_1.4                GenomicFeatures_1.25.12   
##  [3] bitops_1.0-6               tools_3.4.0               
##  [5] zlibbioc_1.19.0            biomaRt_2.29.0            
##  [7] digest_0.6.9               evaluate_0.9              
##  [9] memoise_1.0.0              RSQLite_1.0.0             
## [11] lattice_0.20-33            Matrix_1.2-6              
## [13] graph_1.51.0               DBI_0.4-1                 
## [15] speedglm_0.3-1             withr_1.0.1               
## [17] stringr_1.0.0              grid_3.4.0                
## [19] Biobase_2.33.0             AnnotationDbi_1.35.3      
## [21] XML_3.98-1.4               RBGL_1.49.1               
## [23] BiocParallel_1.7.2         Rsamtools_1.25.0          
## [25] GenomicAlignments_1.9.0    MASS_7.3-45               
## [27] splines_3.4.0              SummarizedExperiment_1.3.2
## [29] stringi_1.0-1              RCurl_1.95-4.8
```

