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
##                                                           SRR1039508 
## "/usr/local/lib/R/site-library/airway/extdata/SRR1039508_subset.bam"
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
## Import genomic features from the file as a GRanges object ... OK
## Prepare the 'metadata' data frame ... OK
## Make the TxDb object ... OK
```

```r
geneids <- keys(txdb, "GENEID")
txdf <- select(txdb, columns=c("TXNAME","TXID"), keys=geneids, keytype="GENEID")
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
We demonstrate the functions nevertheless and plan to create a small
demonstration dataset in the meantime. Robust fitting of these bias
parameters requires 50 or more medium to highly expressed genes.


```r
library(alpine)
library(BSgenome.Hsapiens.UCSC.hg19)
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
##                                                               SRR1039508 
## "/usr/local/lib/R/site-library/airway/extdata/SRR1039508_subset.bam.bai"
```

```r
# here, we can include any combination or functions of bias terms
models <- list("GC"=list(formula="count~ns(gc,knots=gc.knots,Boundary.knots=gc.bk) + gene",
                 offset=c("fraglen")))
# fits one sample at a time
fitpar <- fitModelOverGenes(ebt, bamfiles[1], fragtypes, genome=Hsapiens,
                            models=models,
                            readlength=63, minsize=100, maxsize=300)
fitpar <- list(fitpar) # typically fitpar is a list over samples
names(fitpar) <- names(bamfiles)[1]
```

Visually exploring the bias parameters. These are not robustly fit
in this case because the paucity of reads and genes in the example dataset.


```r
plot(fitpar[[1]]$fraglen.density)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

```r
plotGC(fitpar, m="GC")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-2.png) 

```r
fitpar[[1]]$summary
```

```
## $GC
##                                                     Estimate Std. Error
## (Intercept)                                       -0.3432635  2.1904117
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)1 -0.9413278  2.1268374
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)2 -1.3040727  1.3588469
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)3  1.1721807  4.1833391
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)4  5.5232603  1.6309254
## gene2                                             -1.5840546  0.1182760
## gene3                                             -0.1529903  0.1689949
## gene4                                              0.0987878  0.0969587
## gene5                                             -0.5230994  0.1007034
##                                                    z value Pr(>|z|)
## (Intercept)                                        -0.1567 8.75e-01
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)1  -0.4426 6.58e-01
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)2  -0.9597 3.37e-01
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)3   0.2802 7.79e-01
## ns(gc, knots = gc.knots, Boundary.knots = gc.bk)4   3.3866 7.08e-04
## gene2                                             -13.3929 6.66e-41
## gene3                                              -0.9053 3.65e-01
## gene4                                               1.0189 3.08e-01
## gene5                                              -5.1945 2.05e-07
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
                       subset=TRUE, zerotopos=20, niter=100,
                       lib.sizes=lib.sizes,
                       minsize=100, maxsize=300)
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
## 5.052264e+01 7.392390e+00 1.702922e-31 1.081722e+00 
## 
## $SRR1039508$null$lambda
## 22 23 24 25 
##  1  1  1  1 
## 
## 
## $SRR1039508$GC
## $SRR1039508$GC$theta
##           22           23           24           25 
## 2.851893e+04 4.273922e+03 2.127693e-27 6.299669e+02 
## 
## $SRR1039508$GC$lambda
##          22          23          24          25 
## 0.001770908 0.001751904 0.001685804 0.001720982 
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
## 22 1.278600e+04
## 23 1.916144e+03
## 24 9.539170e-28
## 25 2.824355e+02
## 1  1.625548e-05
## 2  1.783360e+05
## 3  1.029339e+03
## 4  3.041709e+03
## 5  4.853580e+03
## 6  1.732480e+04
## 7  1.737253e+04
## 8  8.496883e-58
## 9  4.497779e+03
## 10 1.679633e+04
## 57 2.880307e+03
## 58 3.894259e+04
## 59 1.209384e+03
## 60 3.816513e-30
## 61 1.251652e-09
## 62 1.064037e+03
## 63 5.079015e-71
```


```r
sessionInfo()
```

```
## R version 3.2.2 (2015-08-14)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 15.10
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
##  [1] knitr_1.11                        BSgenome.Hsapiens.UCSC.hg19_1.4.0
##  [3] BSgenome_1.36.3                   rtracklayer_1.28.10              
##  [5] alpine_0.1.1                      devtools_1.9.1                   
##  [7] GenomicFeatures_1.20.6            AnnotationDbi_1.30.1             
##  [9] Biobase_2.28.0                    GenomicAlignments_1.4.2          
## [11] Rsamtools_1.20.5                  Biostrings_2.36.4                
## [13] XVector_0.8.0                     airway_1.0.0                     
## [15] GenomicRanges_1.20.8              GenomeInfoDb_1.4.3               
## [17] IRanges_2.2.9                     S4Vectors_0.6.6                  
## [19] BiocGenerics_0.14.0              
## 
## loaded via a namespace (and not attached):
##  [1] formatR_1.2.1        compiler_3.2.2       git2r_0.11.0        
##  [4] futile.logger_1.4.1  bitops_1.0-6         futile.options_1.0.0
##  [7] tools_3.2.2          zlibbioc_1.14.0      biomaRt_2.24.1      
## [10] digest_0.6.8         evaluate_0.8         RSQLite_1.0.0       
## [13] memoise_0.2.1        lattice_0.20-33      Matrix_1.2-2        
## [16] DBI_0.3.1            curl_0.9.3           speedglm_0.3-1      
## [19] httr_1.0.0           stringr_1.0.0        grid_3.2.2          
## [22] R6_2.1.1             XML_3.98-1.3         BiocParallel_1.2.22 
## [25] lambda.r_1.1.7       magrittr_1.5         splines_3.2.2       
## [28] MASS_7.3-43          stringi_1.0-1        RCurl_1.95-4.7
```

.
