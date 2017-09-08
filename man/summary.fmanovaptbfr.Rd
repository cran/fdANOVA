\name{summary.fmanovaptbfr}
\alias{summary.fmanovaptbfr}
\title{
Print "fmanovaptbfr" object
}
\description{
Prints the summary of the permutation tests based on a basis function representation for multivariate analysis of variance for functional data, i.e., the W, LH, P and R tests.
}
\usage{
\method{summary}{fmanovaptbfr}(object, \dots)
}
\arguments{
\item{object}{an "\code{fmanovaptbfr}" object.}
\item{...}{additional arguments not used.}
}
\details{
The function prints out the information about the data, the values of test statistics, the p-values and used parameters for the permutation tests based on a basis function representation for FMANOVA problem performed by the \code{\link{fmanova.ptbfr}} function.
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
\seealso{
\code{\link{fmanova.ptbfr}}, \code{\link{print.fmanovaptbfr}}
}
\examples{
# gait data (both features)
library(fda)
gait.data.frame <- as.data.frame(gait)
x.gait <- vector("list", 2)
x.gait[[1]] <- as.matrix(gait.data.frame[, 1:39])
x.gait[[2]] <- as.matrix(gait.data.frame[, 40:78])

# vector of group labels
group.label.gait <- rep(1:3, each = 13)

# the tests based on a basis function representation with default parameters
set.seed(123)
fmanova1 <- fmanova.ptbfr(x.gait, group.label.gait)
summary(fmanova1)
}
\keyword{Printing}
