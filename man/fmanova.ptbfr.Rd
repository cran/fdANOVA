\name{fmanova.ptbfr}
\alias{fmanova.ptbfr}
\title{
Permutation Tests Based on a Basis Function Representation for FMANOVA Problem
}
\description{
Performs the permutation tests based on a basis function representation for multivariate analysis of variance for functional data, i.e., the W, LH, P and R tests. For details of the tests, see Section 2.2 of the vignette file (\code{vignette("fdANOVA", package = "fdANOVA")}).

We consider independent vectors of random functions \deqn{\mathbf{X}_{ij}(t)=(X_{ij1}(t),\dots,X_{ijp}(t))^{\top}\in SP_p(\boldsymbol{\mu}_i,\mathbf{\Gamma}),} \eqn{i=1,\dots,l,} \eqn{j=1,\dots,n_{i}} defined over the interval \eqn{I=[a,b]}, where \eqn{SP_{p}(\boldsymbol{\mu},\boldsymbol{\Gamma})} is a set of \eqn{p}-dimensional stochastic processes with mean vector \eqn{\boldsymbol{\mu}(t)}, \eqn{t\in I} and covariance function \eqn{\boldsymbol{\Gamma}(s, t)}, \eqn{s,t\in I}. Let \eqn{n=n_1+\dots+n_l}. In the multivariate analysis of variance problem for functional data (FMANOVA), we have to test the null hypothesis as follows:
\deqn{
H_0:\boldsymbol{\mu}_1(t)=\dots=\boldsymbol{\mu}_l(t),\ t\in I.
}
The alternative is the negation of the null hypothesis. We assume that each functional observation is observed on a common grid of \eqn{\mathcal{T}} design time points equally spaced in \eqn{I} (see Section 3.1 of the vignette file, \code{vignette("fdANOVA", package = "fdANOVA")}).
}
\usage{
fmanova.ptbfr(x = NULL, group.label, int, B = 1000,
              parallel = FALSE, nslaves = NULL,
              basis = c("Fourier", "b-spline", "own"),
              own.basis, own.cross.prod.mat,
              criterion = c("BIC", "eBIC", "AIC", "AICc", "NO"),
              commonK = c("mode", "min", "max", "mean"),
              minK = NULL, maxK = NULL, norder = 4, gamma.eBIC = 0.5)
}
\arguments{
\item{x}{a list of \eqn{\mathcal{T}\times n} matrices of data, whose each column is a discretized version of a function  and rows correspond to design time points. The \eqn{m}th element of this list contains the data of \eqn{m}th feature, \eqn{m=1,\dots,p}. Its default values is \code{NULL}, because a basis representation of the data can be given instead of raw observations (see the parameter \code{own.basis} below).
}
\item{group.label}{a vector containing group labels.
}
\item{int}{a vector of two elements representing the interval \eqn{I=[a,b]}. When it is not specified, it is determined by a number of design time points.
}
\item{B}{a number of permutation replicates.
}
\item{parallel}{a logical indicating whether to use parallelization.
}
\item{nslaves}{if \code{parallel = TRUE}, a number of slaves. Its default value means that it will be equal to a number of logical processes of a computer used.
}
\item{basis}{a choice of basis of functions used in the basis function representation of the data.
}
\item{own.basis}{if \code{basis = "own"}, a list of length \eqn{p}, whose elements are \eqn{K_m\times n} matrices (\eqn{m=1,\dots,p}) with columns containing the coefficients of the basis function representation of the observations.
}
\item{own.cross.prod.mat}{if \code{basis = "own"}, a \eqn{KM\times KM} cross product matrix corresponding to a basis used to obtain the list \code{own.basis}.
}
\item{criterion}{a choice of information criterion for selecting the optimum value of \eqn{K_m}, \eqn{m=1,\dots,p}. \code{criterion = "NO"} means that \eqn{K_m} are equal to the parameter \code{maxK} defined below. Further remarks about this argument are the same as for the function \code{\link{fanova.tests}}.
}
\item{commonK}{a choice of method for selecting the common value for all observations from the values of \eqn{K_m} corresponding to all processes.
}
\item{minK}{a minimum value of \eqn{K_m}. Further remarks about this argument are the same as for the function \code{\link{fanova.tests}}.
}
\item{maxK}{a maximum value of \eqn{K_m}. Further remarks about this argument are the same as for the function \code{\link{fanova.tests}}.
}
\item{norder}{if \code{basis = "b-spline"}, an integer specifying the order of b-splines.
}
\item{gamma.eBIC}{a \eqn{\gamma\in[0,1]} parameter in the eBIC.
}
}
\value{
A list with class "fmanovaptbfr" containing the following components:
\item{W}{a value of the statistic W.}
\item{pvalueW}{p-value for the W test.}
\item{LH}{a value of the statistic LH.}
\item{pvalueLH}{p-value for the LH test.}
\item{P}{a value of the statistic P.}
\item{pvalueP}{p-value for the P test.}
\item{R}{a value of the statistic R.}
\item{pvalueR}{p-value for the R test,}
the values of parameters used and eventually
\item{data}{a list containing the data given in \code{x}.}
\item{Km}{a vector \eqn{(K_1,\dots,K_p)}.}
\item{KM}{a maximum of a vector \code{Km}.}
}
\references{
Gorecki T, Smaga L (2015). A Comparison of Tests for the One-Way ANOVA Problem for
Functional Data. Computational Statistics 30, 987-1010.

Gorecki T, Smaga L (2017). Multivariate Analysis of Variance for Functional Data. Journal of Applied Statistics 44, 2172-2189.
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fanova.tests}}, \code{\link{fmanova.trp}}, \code{\link{plotFANOVA}}
}
\examples{
# Some of the examples may run some time.

# gait data (both features)
library(fda)
gait.data.frame <- as.data.frame(gait)
x.gait <- vector("list", 2)
x.gait[[1]] <- as.matrix(gait.data.frame[, 1:39])
x.gait[[2]] <- as.matrix(gait.data.frame[, 40:78])

# vector of group labels
group.label.gait <- rep(1:3, each = 13)
\donttest{
# the tests based on a basis function representation with default parameters
set.seed(123)
(fmanova1 <- fmanova.ptbfr(x.gait, group.label.gait))
summary(fmanova1)

# the tests based on a basis function representation with non-default parameters
set.seed(123)
fmanova2 <- fmanova.ptbfr(x.gait, group.label.gait, int = c(0.025, 0.975), B = 5000,
                          basis = "b-spline", criterion = "eBIC", commonK = "mean",
                          minK = 5, maxK = 20, norder = 4, gamma.eBIC = 0.7)
summary(fmanova2)

# the tests based on a basis function representation
# with predefined basis function representation
library(fda)
fbasis <- create.fourier.basis(c(0, nrow(x.gait[[1]])), 17)
own.basis <- vector("list", 2)
own.basis[[1]] <- Data2fd(1:nrow(x.gait[[1]]), x.gait[[1]], fbasis)$coefs
own.basis[[2]] <- Data2fd(1:nrow(x.gait[[2]]), x.gait[[2]], fbasis)$coefs
own.cross.prod.mat <- diag(rep(1, 17))
set.seed(123)
fmanova3 <- fmanova.ptbfr(group.label = group.label.gait,
                          B = 1000, basis = "own",
                          own.basis = own.basis,
                          own.cross.prod.mat = own.cross.prod.mat)
summary(fmanova3)

library(fda)
fbasis <- create.bspline.basis(c(0, nrow(x.gait[[1]])), 20, norder = 4)
own.basis <- vector("list", 2)
own.basis[[1]] <- Data2fd(1:nrow(x.gait[[1]]), x.gait[[1]], fbasis)$coefs
own.basis[[2]] <- Data2fd(1:nrow(x.gait[[2]]), x.gait[[2]], fbasis)$coefs
own.cross.prod.mat <- inprod(fbasis, fbasis)
set.seed(123)
fmanova4 <- fmanova.ptbfr(group.label = group.label.gait,
                          B = 1000, basis = "own",
                          own.basis = own.basis,
                          own.cross.prod.mat = own.cross.prod.mat)
summary(fmanova4)}
\dontshow{
# the tests based on a basis function representation with non-default parameters
fmanova2 <- fmanova.ptbfr(x.gait, group.label.gait, int = c(0.025, 0.975), B = 5,
                          basis = "b-spline", criterion = "eBIC", commonK = "mean",
                          minK = 5, maxK = 20, norder = 4, gamma.eBIC = 0.7)

# the tests based on a basis function representation
# with predefined basis function representation
fbasis <- create.fourier.basis(c(0, nrow(x.gait[[1]])), 17)
own.basis <- vector("list", 2)
own.basis[[1]] <- Data2fd(1:nrow(x.gait[[1]]), x.gait[[1]], fbasis)$coefs
own.basis[[2]] <- Data2fd(1:nrow(x.gait[[2]]), x.gait[[2]], fbasis)$coefs
own.cross.prod.mat <- diag(rep(1, 17))
fmanova3 <- fmanova.ptbfr(group.label = group.label.gait,
                          B = 10, basis = "own",
                          own.basis = own.basis,
                          own.cross.prod.mat = own.cross.prod.mat)

fbasis <- create.bspline.basis(c(0, nrow(x.gait[[1]])), 20, norder = 4)
own.basis <- vector("list", 2)
own.basis[[1]] <- Data2fd(1:nrow(x.gait[[1]]), x.gait[[1]], fbasis)$coefs
own.basis[[2]] <- Data2fd(1:nrow(x.gait[[2]]), x.gait[[2]], fbasis)$coefs
own.cross.prod.mat <- inprod(fbasis, fbasis)
fmanova4 <- fmanova.ptbfr(group.label = group.label.gait,
                          B = 10, basis = "own",
                          own.basis = own.basis,
                          own.cross.prod.mat = own.cross.prod.mat)}
}
\keyword{Multivariate}
\keyword{Test}
