\name{summary.fanovatests}
\alias{summary.fanovatests}
\title{
Print "fanovatests" object
}
\description{
Prints the summary of the FANOVA tests.
}
\usage{
\method{summary}{fanovatests}(object, \dots)
}
\arguments{
\item{object}{an "\code{fanovatests}" object.}
\item{...}{additional arguments not used.}
}
\details{
The function prints out the information about the data, the values of test statistics, the p-values and used parameters for the FANOVA tests performed by the \code{\link{fanova.tests}} function.
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
\seealso{
\code{\link{fanova.tests}}, \code{\link{print.fanovatests}}
}
\examples{
# This example may run some time.

# gait data (the first feature)
library(fda)
gait.data.frame <- as.data.frame(gait)
x.gait <- as.matrix(gait.data.frame[, 1:39])

# vector of group labels
group.label.gait <- rep(1:3, each = 13)
\donttest{
# all FANOVA tests with default parameters
set.seed(123)
fanova1 <- fanova.tests(x.gait, group.label.gait)
summary(fanova1)}
\dontshow{
# only three tests with non-default parameters
fanova2 <- fanova.tests(x.gait, group.label.gait,
                        test = c("FP", "GPF", "Fmaxb"),
                        params = list(paramFP = list(int = c(0.025, 0.975),
                                                     B.FP = 3, basis = "b-spline",
                                                     criterion = "eBIC",
                                                     method = "mean",
                                                     minK = 5, maxK = 20,
                                                     norder = 4, gamma.eBIC = 0.7),
                                      paramFmaxb = 3))
summary(fanova2)}
}
\keyword{Printing}
