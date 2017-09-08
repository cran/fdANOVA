\name{print.fmanovatrp}
\alias{print.fmanovatrp}
\title{
Print "fmanovatrp" object
}
\description{
Prints the brief summary of the testing procedures based on random projections for multivariate analysis of variance for functional data, i.e., the Wp, LHp, Pp and Rp tests.
}
\usage{
\method{print}{fmanovatrp}(x, \dots)
}
\arguments{
\item{x}{an "\code{fmanovatrp}" object.}
\item{...}{additional arguments not used.}
}
\details{
The function prints out the p-values of the tests based on random projections for FMANOVA problem performed by the \code{\link{fmanova.trp}} function.
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
\seealso{
\code{\link{fmanova.trp}}, \code{\link{summary.fmanovatrp}}
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

# the tests based on random projections with the Gaussian white noise generated for projections
set.seed(123)
fmanova1 <- fmanova.trp(x.gait, group.label.gait, k = c(1, 5, 10, 15, 20))
print(fmanova1)
}
\keyword{Printing}
