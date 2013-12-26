\name{plot.LPS}
\alias{plot.LPS}
\title{
  Plot method for LPS objects
}

\description{
  This function plots the distributions of the LPS scores in each group for a fitted \code{\link{LPS}} object.
}

\usage{
  \method{plot}{LPS}(x, y, method=c("Wright", "Radmacher", "exact"), xlim, yaxt = "n", xlab = "LPS",
    ylab, ...)
}

\arguments{
  \item{x}{
    An object of class \code{"LPS"}, as returned by \code{\link{LPS}}.
}
  \item{y}{
    Single character value defining y axis : "density" or (bayesian) "probability".
}
  \item{method}{
    Single character value, the method to use for predictions. See \code{\link{predict.LPS}}.
}
  \item{xlim}{
    To be passed to \code{\link{plot}}, see \code{\link{plot.default}}.
}
  \item{yaxt}{
    To be passed to \code{\link{plot}}, see \code{\link{par}}.
}
  \item{xlab}{
    To be passed to \code{\link{plot}}, see \code{\link{plot.default}}.
}
  \item{ylab}{
    To be passed to \code{\link{plot}}, see \code{\link{plot.default}}.
}
  \item{\dots}{
    Further arguments to be passed to \code{\link{plot}} or \code{\link{par}}.
}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{LPS}}
}

\examples{
  # Data with features in columns
  data(rosenwald)
  group <- rosenwald.cli$group
  expr <- t(rosenwald.expr)
  
  # NA imputation (feature's mean to minimize impact)
  f <- function(x) { x[ is.na(x) ] <- round(mean(x, na.rm=TRUE), 3); x }
  expr <- apply(expr, 2, f)
  
  # Coefficients
  coeff <- LPS.coeff(data=expr, response=group)
  
  # 10 best features model
  m <- LPS(data=expr, coeff=coeff, response=group, k=10)
  
  
  # Distributions of scores in each group
  plot(m, "density")
  
  # Probability for each group along the score axis
  plot(m, "probability", yaxt="s")
}

