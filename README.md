## alpine has moved to Bioconductor (8/18/2016)

Furthermore, the version here is already out-of-date. The version at 
Bioconductor has a simplified interface.

I will still use the GitHub Issues tracker for new feature ideas &
TODOs, but the new website for alpine is:

<http://bioconductor.org/packages/alpine>

New feature requests can be submitted as [Issues](https://github.com/mikelove/alpine/issues) here on GitHub.

Bug reports and any other support questions should be posted to:

<http://support.bioconductor.org>

![alpine](http://mikelove.nfshost.com/img/alpine.jpg)

(the [Sassolungo](https://en.wikipedia.org/wiki/Langkofel) mountain in the Dolomites)

*alpine* is an R/Bioconductor package for modeling and correcting fragment
sequence bias for RNA-seq transcript abundance estimation. 

An example workflow can be found in `vignettes/alpine.Rmd`, or by typing:

```{r}
vignette("alpine")
```

*alpine* is currently designed for un-stranded paired-end RNA-seq data.

A manuscript explaining the *alpine* methods and the background behind fragment 
sequence bias is posted to [bioRxiv](http://biorxiv.org/content/early/2015/08/28/025767).
