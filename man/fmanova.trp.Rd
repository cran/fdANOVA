\name{fmanova.trp}
\alias{fmanova.trp}
\title{
Tests Based on Random Projections for FMANOVA Problem
}
\description{
Performs the testing procedures based on random projections for multivariate analysis of variance for functional data, i.e., the Wp, LHp, Pp and Rp tests. For details of the tests, see Section 2.2 of the vignette file (\code{vignette("fdANOVA", package = "fdANOVA")}).

We consider independent vectors of random functions \deqn{\mathbf{X}_{ij}(t)=(X_{ij1}(t),\dots,X_{ijp}(t))^{\top}\in SP_p(\boldsymbol{\mu}_i,\mathbf{\Gamma}),} \eqn{i=1,\dots,l,} \eqn{j=1,\dots,n_{i}} defined over the interval \eqn{I=[a,b]}, where \eqn{SP_{p}(\boldsymbol{\mu},\boldsymbol{\Gamma})} is a set of \eqn{p}-dimensional stochastic processes with mean vector \eqn{\boldsymbol{\mu}(t)}, \eqn{t\in I} and covariance function \eqn{\boldsymbol{\Gamma}(s, t)}, \eqn{s,t\in I}. Let \eqn{n=n_1+\dots+n_l}. In the multivariate analysis of variance problem for functional data (FMANOVA), we have to test the null hypothesis as follows:
\deqn{
H_0:\boldsymbol{\mu}_1(t)=\dots=\boldsymbol{\mu}_l(t),\ t\in I.
}
The alternative is the negation of the null hypothesis. We assume that each functional observation is observed on a common grid of \eqn{\mathcal{T}} design time points equally spaced in \eqn{I} (see Section 3.1 of the vignette file, \code{vignette("fdANOVA", package = "fdANOVA")}).
}
\usage{
fmanova.trp(x, group.label, k = 30, projection = c("GAUSS", "BM"),
            permutation = FALSE, B = 1000,
            independent.projection.tests = TRUE,
            parallel = FALSE, nslaves = NULL)
}
\arguments{
\item{x}{a list of \eqn{\mathcal{T}\times n} matrices of data, whose each column is a discretized version of a function  and rows correspond to design time points. The \eqn{m}th element of this list contains the data of \eqn{m}th feature, \eqn{m=1,\dots,p}.
}
\item{group.label}{a vector containing group labels.
}
\item{k}{a vector of numbers of projections.
}
\item{projection}{a method of generating Gaussian processes in step 1 of the test based on random projections presented in Section 2 of the vignette file. If \code{projection = "GAUSS"}, the Gaussian white noise is generated as in the function \code{anova.RPm} from the R package \pkg{fda.usc}. In the second case, the Brownian motion is generated.
}
\item{permutation}{a logical indicating whether to compute p-values by permutation method.
}
\item{B}{a number of permutation replicates.
}
\item{independent.projection.tests}{a logical indicating whether to generate the random projections independently or dependently for different elements of vector \code{k}. In the first case, the random projections for each element of vector \code{k} are generated separately, while in the second one, they are generated as chained subsets, e.g., for \code{k = c(5, 10)}, the first 5 projections are a subset of the last 10. The second way of generating random projections is faster than the first one.
}
\item{parallel}{a logical indicating whether to use parallelization.
}
\item{nslaves}{if \code{parallel = TRUE}, a number of slaves. Its default value means that it will be equal to a number of logical processes of a computer used.
}
}
\value{
A list with class "fmanovatrp" containing the following components (|\code{k}| denotes the length of vector \code{k}):
\item{pvalues}{a \eqn{4\times} |\code{k}| matrix of p-values of the tests.}
\item{data.projections}{if \code{independent.projection.tests = TRUE}, a list of length |\code{k}|, whose elements are lists of \eqn{n\times p} matrices of projections of the observations, while when \code{independent.projection.tests = FALSE}, a list of length \eqn{\max}(\code{k}), whose elements are \eqn{n\times p} matrices of projections of the observations,}
and the values of used parameters.
}
\references{
Cuesta-Albertos JA, Febrero-Bande M (2010). A Simple Multiway ANOVA for Functional
Data. Test 19, 537-557.

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
\code{\link{fanova.tests}}, \code{\link{fmanova.ptbfr}}, \code{\link{plotFANOVA}}, \code{\link{plot.fmanovatrp}}
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
# the tests based on random projections with the Gaussian white noise generated for projections
set.seed(123)
fmanova1 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20))
summary(fmanova1)
fmanova2 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20),
                        permutation = TRUE, B = 1000,
                        parallel = TRUE, nslaves = 2)
summary(fmanova2)

# the tests based on random projections with the Brownian motion generated for projections
set.seed(123)
fmanova3 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20), projection = "BM")
summary(fmanova3)
fmanova4 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20), projection = "BM",
                        permutation = TRUE, B = 1000,
                        parallel = TRUE, nslaves = 2)
summary(fmanova4)}
\dontshow{
# the tests based on random projections with the Gaussian white noise generated for projections
fmanova1 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 3))

# the tests based on random projections with the Brownian motion generated for projections
fmanova3 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 3), projection = "BM")}
}
\keyword{Multivariate}
\keyword{Test}
