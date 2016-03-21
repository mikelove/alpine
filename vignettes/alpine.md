<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{alpine}
-->

# Modeling and correcting fragment sequence bias

Here we show a brief example of using the *alpine* package to model
bias parameters and then using those parameters to estimate transcript
abundance. The core *alpine* functions will soon be be wrapped into
convenience functions. First, load some small, subset BAM files from
the *airway* package. 




```r
library(airway)
library(GenomicAlignments)
library(GenomicFeatures)
```


```r
dir <- system.file("extdata", package="airway")
list.files(dir)
```

```
##  [1] "GSE52778_series_matrix.txt"       
##  [2] "Homo_sapiens.GRCh37.75_subset.gtf"
##  [3] "sample_table.csv"                 
##  [4] "SraRunInfo_SRP033351.csv"         
##  [5] "SRR1039508_subset.bam"            
##  [6] "SRR1039508_subset.bam.bai"        
##  [7] "SRR1039509_subset.bam"            
##  [8] "SRR1039512_subset.bam"            
##  [9] "SRR1039513_subset.bam"            
## [10] "SRR1039516_subset.bam"            
## [11] "SRR1039517_subset.bam"            
## [12] "SRR1039520_subset.bam"            
## [13] "SRR1039521_subset.bam"
```

```r
bamfiles <- list.files(dir, "bam$", full=TRUE)
names(bamfiles) <- sub("(.*)_subset.bam","\\1",basename(bamfiles))
bamfiles[1]
```

```
##                                                                                            SRR1039508 
## "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/airway/extdata/SRR1039508_subset.bam"
```

These are reads from a small region.


```r
ga <- readGAlignments(bamfiles[1])
range(ranges(ga))
```

```
## IRanges object with 1 range and 0 metadata columns:
##           start       end     width
##       <integer> <integer> <integer>
##   [1]  11053773  11386194    332422
```

To fit the bias model, we need to identify transcripts which belong to
single-isoform genes.


```r
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
```

```
## Import genomic features from the file as a GRanges object ...
```

```
## OK
```

```
## Prepare the 'metadata' data frame ...
```

```
## OK
```

```
## Make the TxDb object ...
```

```
## OK
```

```r
geneids <- keys(txdb, "GENEID")
txdf <- select(txdb, columns=c("TXNAME","TXID"), keys=geneids, keytype="GENEID")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
# normally, we would pick a set of single isoform genes
tab <- table(txdf$GENEID)
single.tx.genes <- names(tab)[tab == 1]
single.txs <- sort(txdf$TXID[txdf$GENEID %in% single.tx.genes])
# this dataset is too small, so we pick one tx per gene (not recommended).
# we will soon update the vignette to use better example data
txs <- sort(sapply(split(txdf$TXID, txdf$GENEID), `[`, 1))
ebt <- exonsBy(txdb, "tx")
ebt <- ebt[txs]
```

These transcripts should have medium to high coverage.


```r
so <- summarizeOverlaps(ebt, bamfiles, singleEnd=FALSE)
ebt <- ebt[apply(assay(so),1,min) > 50]
so <- summarizeOverlaps(ebt, bamfiles, singleEnd=FALSE)
assay(so)
```

```
##    SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517
## 1         353        271        461        202        439        501
## 12        131        150        180        187        192        261
## 22        113        118        123        120        124        241
## 39        189        229        242        225        234        407
## 45        598        696        743        736        642        983
## 57        346        440        387        591        494        742
##    SRR1039520 SRR1039521
## 1         283        279
## 12        122        247
## 22        103        198
## 39        147        312
## 45        516       1103
## 57        276        613
```

An example of fitting the bias model. Here, we don't have enough data
to properly fit the model because these BAM files are too small of a subset.
We demonstrate the functions and plan to create a small
demonstration dataset in the meantime. Robust fitting of these bias
parameters is best with ~100 medium to highly expressed genes.


```r
library(alpine)
library(BSgenome.Hsapiens.UCSC.hg19)
```

```
## Loading required package: BSgenome
```

```
## Loading required package: rtracklayer
```

```r
seqlevelsStyle(Hsapiens) <- "NCBI" # because these BAMs are NCBI-style
genenames <- names(ebt)
names(genenames) <- genenames
# list of fragment types for each single-isoform gene
fragtypes <- lapply(genenames, function(gene) {
               buildFragtypes(ebt[[gene]], genome=Hsapiens,
               readlength=63, minsize=100, maxsize=300)
             })
indexBam(bamfiles[1])
```

```
##                                                                                                SRR1039508 
## "/Library/Frameworks/R.framework/Versions/3.3/Resources/library/airway/extdata/SRR1039508_subset.bam.bai"
```

The definition of bias models is extremely flexible in *alpine*.
These are defined in a list structure, where each element is a list
with elements `formula` and `offset`. `offset` can be set to `NULL`.
Any kind of R formula is allowed here, which uses information stored
in the elements contained in `fragtypes`. 
The allowable offsets are `fraglen` and/or `vlmm` which
are listed as a character vector. Here we fit a bias model
using fragment length, fragment GC content, and a term for
differences in expression across the genes (`+ gene`).


```r
models <- list("GC"=list(formula="count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + gene",
                         offset=c("fraglen")))
```

Below is an example of full model which is used in the *alpine* paper, but
which we do not fit here. The knots and boundary knots for GC content
and relative position splines are currently fixed internally, and should be
referred to using the following calls:


```r
  "all" = list(formula = "count ~ ns(gc,knots=gc.knots,Boundary.knots=gc.bk) +
  ns(relpos,knots=relpos.knots,Boundary.knots=relpos.bk) +
  GC40.80 + GC40.90 + GC20.80 + GC20.90 +
  gene",
  offset=c("fraglen","vlmm"))
```

The following command then fits bias parameters one sample at a time:


```r
fitpar <- fitModelOverGenes(genes=ebt,
                            bamfile=bamfiles[1],
                            fragtypes=fragtypes,
                            genome=Hsapiens, models=models,
                            readlength=63, minsize=100, maxsize=300)
fitpar <- list(fitpar) # typically fitpar is a list over samples
names(fitpar) <- names(bamfiles)[1]
```

Visually exploring the bias parameters. These are not robustly fit
in this case because the paucity of reads and genes in the example dataset.


```r
plotFragLen(fitpar)
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

```r
plotGC(fitpar, model="GC")
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-2.png)

```r
print(fitpar[[1]]$summary$GC, row.names=FALSE)
```

```
##    Estimate Std. Error  z value Pr(>|z|)
##  -2.0060523  2.1216914  -0.9455 3.44e-01
##   0.7404867  2.0591684   0.3596 7.19e-01
##  -0.2549202  1.3123045  -0.1943 8.46e-01
##   4.6249544  4.0634779   1.1382 2.55e-01
##   5.9316885  1.5189891   3.9050 9.42e-05
##  -1.7058758  0.1157174 -14.7417 3.48e-49
##   0.1088755  0.1543627   0.7053 4.81e-01
##  -0.0186513  0.0936166  -0.1992 8.42e-01
##  -0.6162663  0.0979522  -6.2915 3.14e-10
```

Estimate transcript abundance, first pick a multiple isoform gene.


```r
tab <- table(txdf$GENEID)
mult.tx.genes <- names(tab)[tab > 1]
txs2 <- sort(txdf$TXID[txdf$GENEID %in% mult.tx.genes])
ebt2 <- exonsBy(txdb, "tx")
ebt2 <- ebt2[txs2]
```

For demonstration, pick a gene that has sufficient fragment count.


```r
so <- summarizeOverlaps(ebt2, bamfiles[1])
o <- order(assay(so), decreasing=TRUE)[1:3]
txs2 <- rownames(so)[o]
geneids <- txdf$GENEID[txdf$TXID %in% txs2]
```


```r
models <- list("null"=list(formula=NULL, offset=NULL),
               "GC"=list(formula="count~ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + 0",
                         offset=c("fraglen")))
```


```r
lib.sizes <- 1e6 # any value, this is corrected after counting fragments
names(lib.sizes) <- names(bamfiles)[1]
ebt2 <- exonsBy(txdb, "tx")
# this estimates FPKM for multiple genes and multiple samples
res <- lapply(geneids, function(geneid) {
         txs <- txdf$TXID[txdf$GENEID == geneid]
         ebt <- ebt2[txs]
         estimateTheta(transcripts=ebt, bamfiles=bamfiles[1],
                       fitpar=fitpar, genome=Hsapiens,
                       models=models, readlength=63,
                       minsize=100, maxsize=300,
                       subset=TRUE, zerotopos=20, niter=100,
                       lib.sizes=lib.sizes)
         })
```


```r
res[[1]]
```

```
## $SRR1039508
## $SRR1039508$null
## $SRR1039508$null$theta
##           22           23           24           25 
## 5.051656e+01 7.452425e+00 8.256590e-32 1.080932e+00 
## 
## $SRR1039508$null$lambda
## 22 23 24 25 
##  1  1  1  1 
## 
## 
## $SRR1039508$GC
## $SRR1039508$GC$theta
##           22           23           24           25 
## 4.952758e+03 6.894653e+02 4.265782e-26 1.163922e+02 
## 
## $SRR1039508$GC$lambda
##          22          23          24          25 
## 0.010207977 0.010391124 0.009623262 0.009570531 
## 
## 
## $SRR1039508$count
## [1] 215
```

These estimates are consistent within sample, but need to be scaled
given the total fragment count and the total bias observed over
genes. The `extract` function does this:


```r
mat <- extract(res, model="GC", nsamp=1)
mat
```

```
##      SRR1039508
## 22 1.285411e+04
## 23 1.789399e+03
## 24 1.107117e-25
## 25 3.020776e+02
## 1  8.319942e-06
## 2  1.768663e+05
## 3  1.039145e+03
## 4  3.093056e+03
## 5  5.653784e+03
## 6  1.628469e+04
## 7  1.691459e+04
## 8  9.862229e-59
## 9  4.385137e+03
## 10 1.771738e+04
## 57 2.959845e+03
## 58 3.967781e+04
## 59 1.307408e+03
## 60 2.238736e-30
## 61 1.408310e-09
## 62 1.038066e+03
## 63 2.608835e-73
```


```r
sessionInfo()
```

```
## R version 3.3.0 alpha (2016-03-17 r70348)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
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
##  [1] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.39.4                  
##  [3] rtracklayer_1.31.7                alpine_0.1.3                     
##  [5] GenomicFeatures_1.23.25           AnnotationDbi_1.33.7             
##  [7] GenomicAlignments_1.7.20          Rsamtools_1.23.4                 
##  [9] Biostrings_2.39.12                XVector_0.11.7                   
## [11] airway_0.105.0                    SummarizedExperiment_1.1.22      
## [13] Biobase_2.31.3                    GenomicRanges_1.23.24            
## [15] GenomeInfoDb_1.7.6                IRanges_2.5.40                   
## [17] S4Vectors_0.9.43                  BiocGenerics_0.17.3              
## [19] knitr_1.12.3                      testthat_0.11.0                  
## [21] devtools_1.10.0                   BiocInstaller_1.21.3             
## 
## loaded via a namespace (and not attached):
##  [1] compiler_3.3.0      formatR_1.3         bitops_1.0-6       
##  [4] tools_3.3.0         zlibbioc_1.17.0     biomaRt_2.27.2     
##  [7] digest_0.6.9        evaluate_0.8.3      memoise_1.0.0      
## [10] RSQLite_1.0.0       lattice_0.20-33     Matrix_1.2-4       
## [13] DBI_0.3.1           speedglm_0.3-1      stringr_1.0.0      
## [16] grid_3.3.0          XML_3.98-1.4        BiocParallel_1.5.20
## [19] magrittr_1.5        MASS_7.3-45         splines_3.3.0      
## [22] stringi_1.0-1       RCurl_1.95-4.8      crayon_1.3.1
```

.
