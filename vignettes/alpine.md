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
##                                                      SRR1039508 
## "/home/love/bin/R/library/airway/extdata/SRR1039508_subset.bam"
```

These are reads from a small region.


```r
ga <- readGAlignments(bamfiles[1])
range(ranges(ga))
```

```
## IRanges of length 1
##        start      end  width
## [1] 11053773 11386194 332422
```

To fit the bias model, we need to identify transcripts which belong to
single-isoform genes.


```r
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs=character())
```

```
## Prepare the 'metadata' data frame ... metadata: OK
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
# this dataset is too small, so we pick one tx per gene (not recommended)
txs <- sort(sapply(split(txdf$TXID, txdf$GENEID), `[`, 1))
ebt <- exonsBy(txdb, "tx")
ebt <- ebt[txs]
```

These transcripts should have medium to high coverage.


```r
so <- summarizeOverlaps(ebt, bamfiles, singleEnd=FALSE)
```

 
 

```r
ebt <- ebt[apply(assay(so),1,min) > 50]
so <- summarizeOverlaps(ebt, bamfiles, singleEnd=FALSE)
```

 
 

```r
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
We demonstrate the functions nevertheless and plan to create a small
demonstration dataset in the meantime. Robust fitting of these bias
parameters requires 50 or more medium to highly expressed genes.


```r
library(alpine)
library(BSgenome.Hsapiens.UCSC.hg19)
```

```
## Loading required package: BSgenome
## Loading required package: rtracklayer
```

```r
seqlevelsStyle(Hsapiens) <- "NCBI" # because the BAMs are NCBI-style
genenames <- names(ebt)
names(genenames) <- genenames
# list of fragment types for each single-isoform gene
fragtypes <- lapply(genenames, function(gene) {
               buildFragtypesFromExons(ebt[[gene]], genome=Hsapiens,
               readlength=63, minsize=100, maxsize=300)
             })
indexBam(bamfiles[1])
```

```
##                                                          SRR1039508 
## "/home/love/bin/R/library/airway/extdata/SRR1039508_subset.bam.bai"
```

```r
# here, we can include many kinds of modeling terms
models <- list("GC"=list(formula="count~ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + gene",
                 offset=c("fraglen")))
# fits one sample at a time
fitpar <- fitModelOverGenes(ebt, bamfiles[1], fragtypes, genome=Hsapiens,
                            models=models,
                            readlength=63, minsize=100, maxsize=300)
fitpar <- list(fitpar)
names(fitpar) <- names(bamfiles)[1]
```

Visually exploring the bias parameters. These are not robustly fit
in this case because the paucity of reads and genes in the example dataset.


```r
plot(fitpar[[1]]$fraglen.density)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

```r
plotOrder0(fitpar[[1]]$vlmm.fivep$order0)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-2.png) 

```r
plotGC(fitpar, m="GC")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-3.png) 

```r
unname(fitpar[[1]]$coefs[["GC"]])
```

```
## [1] -3.40302829  2.28981789  0.53499983  7.85456827  7.51507564 -1.79022474
## [7] -0.04054359 -0.13496314 -0.75902261
```

```r
fitpar[[1]]$summary
```

```
## $GC
##                                                     Estimate Std. Error
## (Intercept)                                       -3.4030283  2.1299511
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)1  2.2898179  2.0655535
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)2  0.5349998  1.3195946
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)3  7.8545683  4.0772535
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)4  7.5150756  1.5228968
## gene2                                             -1.7902247  0.1131065
## gene3                                             -0.0405436  0.1534809
## gene4                                             -0.1349631  0.0902490
## gene5                                             -0.7590226  0.0956312
##                                                    z value Pr(>|z|)
## (Intercept)                                        -1.5977 1.10e-01
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)1   1.1086 2.68e-01
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)2   0.4054 6.85e-01
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)3   1.9264 5.40e-02
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)4   4.9347 8.03e-07
## gene2                                             -15.8278 2.00e-56
## gene3                                              -0.2642 7.92e-01
## gene4                                              -1.4955 1.35e-01
## gene5                                              -7.9370 2.07e-15
```

Estimate transcript abundance, first pick a multiple isoform gene.


```r
tab <- table(txdf$GENEID)
mult.tx.genes <- names(tab)[tab > 1]
txs <- sort(txdf$TXID[txdf$GENEID %in% mult.tx.genes])
ebt <- exonsBy(txdb, "tx")
ebt <- ebt[txs]
```

For demonstration, pick a gene that has sufficient fragment count.


```r
so <- summarizeOverlaps(ebt, bamfiles[1])
```

 

```r
tx <- rownames(so)[which.max(assay(so))]
geneid <- txdf$GENEID[txdf$TXID == tx]
txs <- txdf$TXID[txdf$GENEID == geneid]
ebt <- exonsBy(txdb, "tx")
ebt <- ebt[txs]
```


```r
models <- list("null"=list(formula=NULL, offset=NULL),
               "GC"=list(formula="count~ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + 0",
                 offset=c("fraglen")))
```


```r
lib.sizes <- 30e6 # must be pre-calculated
names(lib.sizes) <- names(bamfiles)[1]
# for efficiency, should be run with all models and all samples at once
res <- estimateTheta(transcripts=ebt, bamfiles=bamfiles[1],
                     fitpar=fitpar, genome=Hsapiens,
                     models=models, readlength=63,
                     subset=TRUE, zerotopos=20, niter=100,
                     lib.sizes=lib.sizes,
                     minsize=100, maxsize=300)
```

These estimates are consistent within sample, but should be scaled
to the null model, and then calibrated across sample using the
median-ratio method. (Functions to come).


```r
res
```

```
## $SRR1039508
## $SRR1039508$null
##           57           58           59           60           61 
## 4.297594e-01 4.889635e+00 2.461233e-01 2.921742e-31 2.142875e-08 
##           62           63 
## 1.526115e-01 1.653457e-68 
## 
## $SRR1039508$GC
##           57           58           59           60           61 
## 7.715833e+00 1.131809e+02 3.826664e+00 2.361135e-34 4.359618e-12 
##           62           63 
## 3.141146e+00 9.203357e-76
```


```r
sessionInfo()
```

```
## R Under development (unstable) (2015-06-21 r68565)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: Ubuntu 15.04
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
##  [1] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.37.3                  
##  [3] rtracklayer_1.29.12               alpine_0.1                       
##  [5] GenomicFeatures_1.21.13           AnnotationDbi_1.31.17            
##  [7] GenomicAlignments_1.5.11          Rsamtools_1.21.13                
##  [9] Biostrings_2.37.2                 XVector_0.9.1                    
## [11] airway_0.103.1                    SummarizedExperiment_0.3.2       
## [13] Biobase_2.29.1                    GenomicRanges_1.21.16            
## [15] GenomeInfoDb_1.5.8                IRanges_2.3.14                   
## [17] S4Vectors_0.7.10                  BiocGenerics_0.15.3              
## [19] testthat_0.10.0                   devtools_1.8.0                   
## [21] knitr_1.10.5                      BiocInstaller_1.19.14            
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.11.6          formatR_1.2          git2r_0.10.1        
##  [4] futile.logger_1.4.1  bitops_1.0-6         futile.options_1.0.0
##  [7] tools_3.3.0          zlibbioc_1.15.0      biomaRt_2.25.1      
## [10] digest_0.6.8         lattice_0.20-31      evaluate_0.7        
## [13] memoise_0.2.1        RSQLite_1.0.0        Matrix_1.2-2        
## [16] DBI_0.3.1            curl_0.9             speedglm_0.3        
## [19] stringr_1.0.0        xml2_0.1.1           rversions_1.0.1     
## [22] grid_3.3.0           XML_3.98-1.3         BiocParallel_1.3.34 
## [25] lambda.r_1.1.7       magrittr_1.5         splines_3.3.0       
## [28] MASS_7.3-42          stringi_0.5-5        RCurl_1.95-4.7      
## [31] crayon_1.3.0
```

.
