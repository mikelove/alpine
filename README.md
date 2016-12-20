# alpine

![alpine](http://mikelove.nfshost.com/img/alpine.jpg)

(the [Sassolungo](https://en.wikipedia.org/wiki/Langkofel) mountain in the Dolomites)

*alpine* is an R/Bioconductor package for modeling and correcting fragment
sequence bias for RNA-seq transcript abundance estimation. In
particular, it is the first method of its kind to take into account
sample-specific dependence of RNA-seq fragments on their GC content.

An example workflow can be found in `vignettes/alpine.Rmd`, or by typing:

```{r}
vignette("alpine")
```

*alpine* is designed for un-stranded paired-end RNA-seq data.

A paper explaining the *alpine* methods was published in December 2016:

<http://dx.doi.org/10.1038/nbt.3682>
