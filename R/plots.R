#' Plot the fragment GC bias over samples
#' 
#' Plots smooth curves of the log fragment rate over fragment GC content.
#'
#' @param fitpar a list of the output of \link{fitBiasModels} over samples
#' @param model the name of one of the models 
#' @param col a vector of colors
#' @param lty a vector of line types
#' @param ylim the y limits for the plot
#' @param gc.range a numeric of length two,
#' the range of the fragment GC content. By default,
#' [.2,.8] for plotting and [0,1] for returning a matrix
#' @param return.type a numeric, either
#' 0: make a plot,
#' 1: skip the plot and return a matrix of log fragment rate,
#' 2: skip the plot and return a matrix of probabilities 
#'
#' @return Either plot, or if \code{return.type} is 1 or 2, a matrix
#' 
#' @examples
#'
#' # fitpar was fit using identical code
#' # as found in the vignette, except with
#' # 25 genes, and with fragment size in 80-350 bp
#' data(preprocessedData)
#' perf <- rep(1:2, each=2)
#' plotGC(fitpar, "all", col=perf)
#' 
#' @export
plotGC <- function(fitpar, model, col, lty, ylim,
                   gc.range=NULL, return.type=0) {

  if (is.null(gc.range)) {
    gc.range <- if (return.type == 0) {
      c(.2, .8)
    } else {
      c(0,1)
    }
  }
  stopifnot(length(gc.range) == 2)
  stopifnot(length(return.type) == 1 & return.type %in% 0:2)

  # just a single sample?
  if ("models" %in% names(fitpar)) {
    fitpar <- list(fitpar)
  }
  
  knots <- fitpar[[1]][["model.params"]]$gc.knots
  bk <- fitpar[[1]][["model.params"]]$gc.bk
  
  n <- length(knots)
  coef.nms <- names(fitpar[[1]][["coefs"]][[model]])
  coef.idx <- c(grep("\\(Intercept\\)",coef.nms), grep("ns\\(gc", coef.nms))
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[model]][coef.idx])
  gene.coefs <- lapply(fitpar, function(elem)
    elem[["coefs"]][[model]][ grep("gene", names(elem[["coefs"]][[model]])) ])
  # new intercept: the average of the intercept + gene coefficients
  coefmat[1,] <- coefmat[1,] + sapply(gene.coefs, mean)
  z <- seq(from=gc.range[1],to=gc.range[2],length=101)
  x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk))
  logpred <- x %*% coefmat
  rownames(logpred) <- as.character(z)
  if (return.type == 1) {
    return(logpred)
  } else if (return.type == 2) {
    probmat <- exp(logpred)
    probmat <- sweep(probmat, 2, apply(probmat, 2, max), "/")
    return(probmat)
  }
  if (missing(ylim)) {
    ylim <- c(min(logpred),max(logpred))
  }
  plot(0,0,type="n",xlim=c(gc.range[1],gc.range[2]),ylim=ylim,
       ylab="log fragment rate", xlab="fragment GC content",
       main="fragment sequence bias")
  if (missing(col)) {
    col <- rep("black", ncol(logpred))
  }
  if (missing(lty)) {
    lty <- rep(1, ncol(logpred))
  }
  for (i in 1:ncol(logpred)) {
    lines(z, logpred[,i], col=col[i], lwd=2, lty=lty[i])
  }
}

#' Plot relative position bias over samples
#'
#' Plots the smooth curves of log fragment rate over relative position.
#' 
#' @param fitpar a list of the output of \link{fitBiasModels} over samples
#' @param model the name of one of the models 
#' @param col a vector of colors
#' @param lty a vector of line types
#' @param ylim the y limits for the plot
#'
#' @return plot
#' 
#' @examples
#'
#' # fitpar was fit using identical code
#' # as found in the vignette, except with
#' # 25 genes, and with fragment size in 80-350 bp
#' data(preprocessedData)
#' perf <- rep(1:2, each=2)
#' plotRelPos(fitpar, "all", col=perf)
#' 
#' @export
plotRelPos <- function(fitpar, model, col, lty, ylim) {

  # just a single sample?
  if ("models" %in% names(fitpar)) {
    fitpar <- list(fitpar)
  }

  knots <- fitpar[[1]][["model.params"]]$relpos.knots
  bk <- fitpar[[1]][["model.params"]]$relpos.bk

  n <- length(knots)
  coef.nms <- names(fitpar[[1]][["coefs"]][[model]])
  coef.idx <- c(grep("\\(Intercept\\)",coef.nms), grep("ns\\(relpos", coef.nms))
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[model]][coef.idx])
  z <- seq(from=0,to=1,length=101)
  x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk))
  logpred <- x %*% coefmat
  logpred <- scale(logpred, scale=FALSE)
  if (missing(ylim)) {
    ylim <- c(min(logpred),max(logpred))
  }
  plot(0,0,type="n",xlim=c(0,1),ylim=ylim,
       ylab="log fragment rate", xlab="5' -- position in transcript -- 3'",
       main="relative position bias")
  if (missing(col)) {
    col <- rep("black", ncol(logpred))
  }
  if (missing(lty)) {
    lty <- rep(1, ncol(logpred))
  }
  for (i in 1:ncol(logpred)) {
    lines(z, logpred[,i], col=col[i], lwd=2, lty=lty[i])
  }
}

#' Plot fragment length distribution over samples
#'
#' Plots the fragment length distribution.
#' 
#' @param fitpar a list of the output of \link{fitBiasModels} over samples
#' @param col a vector of colors
#' @param lty a vector of line types
#'
#' @return plot
#' 
#' @examples
#' 
#' # fitpar was fit using identical code
#' # as found in the vignette, except with
#' # 25 genes, and with fragment size in 80-350 bp
#' data(preprocessedData)
#' perf <- rep(1:2, each=2)
#' plotFragLen(fitpar, col=perf)
#' 
#' @export
plotFragLen <- function(fitpar, col, lty) {

  # just a single sample?
  if ("models" %in% names(fitpar)) {
    fitpar <- list(fitpar)
  }

  if (missing(col)) {
    col <- rep("black", length(fitpar))
  }
  if (missing(lty)) {
    lty <- rep(1, length(fitpar))
  }
  ymax <- max(sapply(fitpar, function(x) max(x$fraglen.density$y)))
  plot(fitpar[[1]]$fraglen.density, ylim=c(0,ymax*1.1),
       xlab="fragment length", ylab="density", main="fragment length distribution",
       col=col[1], lty=lty[1], lwd=2)
  for (i in seq_along(fitpar)[-1]) {
    lines(fitpar[[i]]$fraglen.density, col=col[i], lty=lty[i], lwd=2)
  }
}

#' Plot parameters of the variable length Markov model (VLMM) for read starts
#'
#' This function plots portions of the Cufflinks VLMM for read start bias.
#' The natural log of observed over expected is shown, such that 0
#' indicates no contribution of a position to the read start bias.
#' As the variable lenght Markov model has different dependencies for different
#' positions (see Roberts et al, 2011), it is difficult
#' to show all the 744 parameters simultaneously. Instead this function
#' offers to show the 0-order terms for all positions, or the 1st and 2nd
#' order terms for selected positions within the read start sequence.
#' For the 1- and 2-order terms, the log bias is shown for each nucleotide
#' (A,C,T,G) given the previous nucleotide (1-order) or di-nucleotide (2-order).
#'
#' @references
#'
#' Roberts et al, "Improving RNA-Seq expression estimates by correcting for fragment bias"
#' Genome Biology (2011) doi:101186/gb-2011-12-3-r22
#' 
#' @param order0 the "order0" element of the list named "vlmm.fivep" or "vlmm.threep"
#' within the list that is the output of \link{fitBiasModels}
#' @param order1 as for "order0" but "order1"
#' @param pos1 the position of the 1st order VLMM to plot
#' @param order2 as for "order0" but "order2"
#' @param pos2 the position of the 2nd order VLMM to plot
#' @param ... parameters passed to \code{plot}
#'
#' @return plot
#' 
#' @examples
#'
#' # fitpar was fit using identical code
#' # as found in the vignette, except with
#' # 25 genes, and with fragment size in 80-350 bp
#' data(preprocessedData)
#' plotOrder0(fitpar[[1]][["vlmm.fivep"]][["order0"]])
#' plotOrder1(fitpar[[1]][["vlmm.fivep"]][["order1"]], pos1=5:19)
#' plotOrder2(fitpar[[1]][["vlmm.fivep"]][["order2"]], pos2=8:17)
#' 
#' @export
plotOrder0 <- function(order0, ...) {
  dna.letters <- c("A","C","G","T")
  mat <- log(order0$obs/order0$expect)
  xpos <- -8:12
  dna.cols <- c("green3","blue3","orange3","red3")
  plot(0,0,xlim=c(-8,12),type="n",xlab="position",ylab="log(observed / expected)", ...)
  for (i in 1:4) {
    points(xpos, mat[i,], col=dna.cols[i], type="b", lwd=2)
  }
  abline(v=0, h=0)
  legend("topright",dna.letters,pch=1,lty=1,col=dna.cols,bg="white")
}

#' @describeIn plotOrder0 Plot first order parameters for a position
#' @export
plotOrder1 <- function(order1, pos1) {
  dna.letters <- c("A","C","G","T")
  order <- 1
  npos <- length(pos1)
  dna.cols <- c("green3","blue3","orange3","red3")
  par(mfrow=c(1,npos+1),mar=c(5,1,3,1))
  for (i in seq_len(npos)) {
    plot(as.vector(log(order1$obs[,,i]/order1$expect)), rev(seq_len(4 * 4^order)),
         col=rep(dna.cols,each=4), xlim=c(-1,1),
         ylab="",xlab="",yaxt="n",main=pos1[i] - 10 + 1,cex=2)
    abline(v=0)
  }
  plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  alpha <- alphafun(dna.letters, order-1)
  legend("center",alpha,pch=1,col=dna.cols,cex=2,title="prev")
}

#' @describeIn plotOrder0 Plot second order parameters for a position
#' @export
plotOrder2 <- function(order2, pos2) {
  dna.letters <- c("A","C","G","T")
  order <- 2
  npos <- length(pos2)
  dna.cols <- c("green3","blue3","orange3","red3")
  par(mfrow=c(1,npos+1),mar=c(5,1,3,1))
  for (i in 1:npos) {
    plot(as.vector(log(order2$obs[,,i]/order2$expect)), rev(seq_len(4 * 4^order)),
         col=rep(dna.cols,each=4),pch=rep(1:4,each=16), xlim=c(-1,1),
         ylab="",xlab="",yaxt="n",main=pos2[i] - 10 + 1,cex=2)
    abline(v=0)
  }
  plot(0,0,type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  alpha <- alphafun(dna.letters, order-1)
  legend("center",alpha,pch=rep(1:4,each=4),col=dna.cols,cex=2,title="prev")
}

#' Simple segments plot for GRangesList
#' 
#' Simple segments plot for GRangesList
#'
#' @param grl GRangesList object
#' @param ... passed to plot
#'
#' @return plot
#' 
#' @examples
#' 
#' library(GenomicRanges)
#' grl <- GRangesList(GRanges("1",IRanges(c(100,200,300),width=50)),
#'                    GRanges("1",IRanges(c(100,300),width=c(75,50))),
#'                    GRanges("1",IRanges(c(100,200,400),width=c(75,50,50))),
#'                    GRanges("1",IRanges(c(200,300,400),width=50)))
#' plotGRL(grl)
#' 
#' @export
plotGRL <- function(grl, ...) {
  df <- as.data.frame(grl)
  plot(0, 0, xlim=range(c(df$start,df$end)), ylim=c(1,max(df$group)),
       type="n", xlab="position", ylab="group", ...)
  segments(df$start, df$group, df$end, df$group, lwd=3)
}
