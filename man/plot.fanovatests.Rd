\name{plot.fanovatests}
\alias{plot.fanovatests}
\title{
Plot P-values of Tests Based on Random Projections for FANOVA Problem
}
\description{
The p-values of the tests based on random projections for the one-way analysis of variance for (univariate) functional data against the number of projections are plotted.
}
\usage{
\method{plot}{fanovatests}(x, y, \dots)
}
\arguments{
\item{x}{an "\code{fanovatests}" object. More precisely, a result of the function \code{\link{fanova.tests}} for the standard tests based on random projections.}
\item{y}{an "\code{fanovatests}" object. More precisely, a result of the function \code{\link{fanova.tests}} for the permutation tests based on random projections.}
\item{...}{additional arguments not used.}
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
\note{
We can use only one of the arguments \code{x} and \code{y}, or both simultaneously.
}
\seealso{
\code{\link{fanova.tests}}, \code{\link{plot.fmanovatrp}}
}
\examples{
# Some of the examples may run some time.

# gait data (the first feature)
library(fda)
gait.data.frame <- as.data.frame(gait)
x.gait <- as.matrix(gait.data.frame[, 1:39])

# vector of group labels
group.label.gait <- rep(1:3, each = 13)
\donttest{
set.seed(123)
fanova4 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30), B.TRP = 1000)))
set.seed(123)
fanova5 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30),
                                                      permutation = TRUE, B.TRP = 1000)))
plot(x = fanova4)
plot(y = fanova5)
plot(x = fanova4, y = fanova5)

set.seed(123)
fanova6 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30), projection = "BM",
                                                      B.TRP = 1000)))
set.seed(123)
fanova7 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30), projection = "BM",
                                                      permutation = TRUE, B.TRP = 1000)))
plot(x = fanova6)
plot(y = fanova7)
plot(x = fanova6, y = fanova7)}
\dontshow{
fanova4 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        params = list(paramTRP = list(k = c(1, 3), B.TRP = 3)))
fanova5 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        params = list(paramTRP = list(k = c(1, 3),
                                                      permutation = TRUE, B.TRP = 3)))
plot(x = fanova4)
plot(y = fanova5)
plot(x = fanova4, y = fanova5)}
}
\keyword{Plot}
