\name{plot.fmanovatrp}
\alias{plot.fmanovatrp}
\title{
Plot P-values of Tests Based on Random Projections for FMANOVA Problem
}
\description{
The p-values of the tests based on random projections for multivariate analysis of variance for functional data against the number of projections are plotted.
}
\usage{
\method{plot}{fmanovatrp}(x, y, withoutRoy = FALSE, \dots)
}
\arguments{
\item{x}{an "\code{fmanovatrp}" object. More precisely, a result of the function \code{\link{fmanova.trp}} for the standard tests based on random projections.}
\item{y}{an "\code{fmanovatrp}" object. More precisely, a result of the function \code{\link{fmanova.trp}} for the permutation tests based on random projections.}
\item{withoutRoy}{a logical indicating whether to plot the p-values of the Rp test.}
\item{...}{additional arguments not used.}
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
\note{
We can use only one of the arguments \code{x} and \code{y}, or both simultaneously.
}
\seealso{
\code{\link{fmanova.trp}}, \code{\link{plot.fanovatests}}
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
fmanova2 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20),
                        permutation = TRUE, B = 1000, parallel = TRUE, nslaves = 2)
plot(x = fmanova1)
plot(x = fmanova1, withoutRoy = TRUE)
plot(y = fmanova2)
plot(x = fmanova1, y = fmanova2)
plot(x = fmanova1, y = fmanova2, withoutRoy = TRUE)

# the tests based on random projections with the Brownian motion generated for projections
set.seed(123)
fmanova3 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20), projection = "BM")
fmanova4 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20), projection = "BM",
                        permutation = TRUE, B = 1000, parallel = TRUE, nslaves = 2)
plot(x = fmanova3)
plot(x = fmanova3, withoutRoy = TRUE)
plot(y = fmanova4)
plot(x = fmanova3, y = fmanova4)
plot(x = fmanova3, y = fmanova4, withoutRoy = TRUE)}
\dontshow{
# the tests based on random projections with the Gaussian white noise generated for projections
fmanova1 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 3))
fmanova2 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 3),
                        permutation = TRUE, B = 5)
plot(x = fmanova1)
plot(y = fmanova2)
plot(x = fmanova1, y = fmanova2)
plot(x = fmanova1, y = fmanova2, withoutRoy = TRUE)

# the tests based on random projections with the Brownian motion generated for projections
fmanova3 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 3), projection = "BM")
fmanova4 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 3), projection = "BM",
                        permutation = TRUE, B = 5)
plot(x = fmanova3)
plot(y = fmanova4)
plot(x = fmanova3, y = fmanova4)
plot(x = fmanova3, y = fmanova4, withoutRoy = TRUE)}
}
\keyword{Plot}
