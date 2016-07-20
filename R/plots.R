#' Plot the fragment GC bias over samples
#' 
#' Plots smooth curves of the log fragment rate over fragment GC content.
#'
#' @param fitpar a list of the output of \link{fitBiasModels} over samples
#' @param model the name of one of the models 
#' @param col a vector of colors
#' @param lty a vector of line types
#' @param ylim the y limits for the plot
#' @param knots the knots for the spline
#' @param bk the boundary knots for the spline
#'
#' @examples
#'
#' # fitpar was fit using identical code
#' # as found in the vignette, except with
#' # 25 genes, and with fragment size in 80-350 bp
#' data(fitpar)
#' perf <- rep(1:2, each=2)
#' plotGC(fitpar, "all", col=perf)
#' 
#' @export
plotGC <- function(fitpar, model, col, lty, ylim, knots=c(.4,.5,.6), bk=c(0,1)) {
  n <- length(knots)
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[model]][1:(n+2)])
  genecoefs <- lapply(fitpar, function(elem)
    elem[["coefs"]][[model]][ grep("gene", names(elem[["coefs"]][[model]])) ])
  # new intercept: the average of the intercept + gene coefficients
  coefmat[1,] <- coefmat[1,] + sapply(genecoefs, mean)
  z <- seq(from=.2,to=.8,length=101)
  x <- model.matrix(~ ns(z, knots=knots, Boundary.knots=bk))
  logpred <- x %*% coefmat
  if (missing(ylim)) {
    ylim <- c(min(logpred),max(logpred))
  }
  plot(0,0,type="n",xlim=c(.2,.8),ylim=ylim,
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
#' @param knots the knots for the spline
#' @param bk the boundary knots for the spline
#' 
#' @export
plotRelPos <- function(fitpar, model, col, lty, ylim, knots=c(.25,.5,.75), bk=c(0,1)) {
  n <- length(knots)
  # assuming same number of knots for GC and for relpos
  coefmat <- sapply(fitpar, function(elem) elem[["coefs"]][[model]][c(1,(3+n):(3+2*n))])
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
#' @export
plotFragLen <- function(fitpar, col, lty) {
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
#' As the variable lenght Markov model has different dependencies for different
#' positions (see Roberts et al, 2011), it is difficult
#' to show all the 744 parameters simultaneously. Instead this function
#' offers to show the 0-order terms for all positions, or the 1st and 2nd
#' order terms for selected positions within the read start sequence.
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
  abline(v=0, h=0, col=rgb(0,0,0,.3))
  legend("topright",dna.letters,pch=1,lty=1,col=dna.cols,bg="white")
}

#' @describeIn plotOrder0 Plot first order parameters for a position
#' @export
plotOrder1 <- function(order1, pos1) {
  order <- 1
  npos <- length(pos1)
  dna.cols <- c("green3","blue3","orange3","red3")
  mypar(1,npos+1,mar=c(5,1,3,1))
  for (i in 1:npos) {
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
  order <- 2
  npos <- length(pos2)
  dna.cols <- c("green3","blue3","orange3","red3")
  mypar(1,npos+1,mar=c(5,1,3,1))
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
#' @export
plotGRL <- function(grl, ...) {
  df <- as.data.frame(grl)
  plot(0, 0, xlim=range(c(df$start,df$end)), ylim=c(1,max(df$group)),
       type="n", xlab="position", ylab="group", ...)
  segments(df$start, df$group, df$end, df$group, lwd=3)
}
