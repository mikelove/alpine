% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predictCoverage}
\alias{predictCoverage}
\title{Predict coverage for a single-isoform gene}
\usage{
predictCoverage(gene, bam.files, fitpar, genome, models, readlength, minsize,
  maxsize)
}
\arguments{
\item{gene}{a GRangesList with the exons of different genes}

\item{bam.files}{a character string pointing to indexed BAM files}

\item{fitpar}{the output of running \code{\link{fitBiasModels}}}

\item{genome}{a BSgenome object}

\item{models}{a list describing the models, see \code{link{fitBiasModels}}}

\item{readlength}{the read length}

\item{minsize}{the minimum fragment length to model}

\item{maxsize}{the maximum fragment length to model}
}
\value{
TODO
}
\description{
Predict coverage for a single-isoform gene given
fitted bias parameters in a set of models,
and compare to the observed fragment coverage.
}
