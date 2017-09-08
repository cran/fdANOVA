\name{plotFANOVA}
\alias{plotFANOVA}
\title{
Plot Univariate Functional Data
}
\description{
Univariate functional observations with or without indication of groups as well as mean functions of samples are plotted. We assume that \eqn{n} univariate functional observations are observed on a common grid of \eqn{\mathcal{T}} design time points equally spaced in \eqn{I=[a,b]} (see Section 3.1 of the vignette file, \code{vignette("fdANOVA", package = "fdANOVA")}).
}
\usage{
plotFANOVA(x, group.label = NULL, int = NULL, separately = FALSE,
           means = FALSE, smooth = FALSE, ...)
}
\arguments{
\item{x}{a \eqn{\mathcal{T}\times n} matrix of data, whose each column is a discretized version of a function and rows correspond to design time points.
}
\item{group.label}{a character vector containing group labels. Its default value means that all functional observations are drawn without division into groups.
}
\item{int}{a vector of two elements representing the interval \eqn{I=[a,b]}. When it is not specified, it is determined by a number of design time points.
}
\item{separately}{a logical indicating how groups are drawn. If \code{separately = FALSE}, groups are drawn on one plot by different colors. When \code{separately = TRUE}, they are depicted in different panels.
}
\item{means}{a logical indicating whether to plot only group mean functions.
}
\item{smooth}{a logical indicating whether to plot reconstructed data via smoothing splines instead of raw data.
}
\item{...}{additional arguments not used.}
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
\seealso{
\code{\link{fanova.tests}}, \code{\link{fmanova.ptbfr}}, \code{\link{fmanova.trp}}
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
plotFANOVA(x = x.gait[[1]], int = c(0.025, 0.975))
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975))
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), separately = TRUE)
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), means = TRUE)

plotFANOVA(x = x.gait[[1]], int = c(0.025, 0.975), smooth = TRUE)
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), smooth = TRUE)
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), separately = TRUE, smooth = TRUE)
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), means = TRUE, smooth = TRUE)

plotFANOVA(x = x.gait[[2]], int = c(0.025, 0.975))
plotFANOVA(x = x.gait[[2]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975))
plotFANOVA(x = x.gait[[2]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), separately = TRUE)
plotFANOVA(x = x.gait[[2]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), means = TRUE)

plotFANOVA(x = x.gait[[2]], int = c(0.025, 0.975), smooth = TRUE)
plotFANOVA(x = x.gait[[2]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), smooth = TRUE)
plotFANOVA(x = x.gait[[2]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), separately = TRUE, smooth = TRUE)
plotFANOVA(x = x.gait[[2]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), means = TRUE, smooth = TRUE)

# Canadian Weather data (both features)
library(fda)
x.CW <- vector("list", 2)
x.CW[[1]] <- CanadianWeather$dailyAv[,,1]
x.CW[[2]] <- CanadianWeather$dailyAv[,,2]

# vector of group labels
group.label.CW <- rep(c("Eastern", "Western", "Northern"), c(15, 15, 5))

plotFANOVA(x = x.CW[[1]])
plotFANOVA(x = x.CW[[1]], group.label = as.character(group.label.CW))
plotFANOVA(x = x.CW[[1]], group.label = as.character(group.label.CW),
           separately = TRUE)
plotFANOVA(x = x.CW[[1]], group.label = as.character(group.label.CW),
           means = TRUE)

plotFANOVA(x = x.CW[[1]], smooth = TRUE)
plotFANOVA(x = x.CW[[1]], group.label = as.character(group.label.CW),
           smooth = TRUE)
plotFANOVA(x = x.CW[[1]], group.label = as.character(group.label.CW),
           separately = TRUE, smooth = TRUE)
plotFANOVA(x = x.CW[[1]], group.label = as.character(group.label.CW),
           means = TRUE, smooth = TRUE)

plotFANOVA(x = x.CW[[2]])
plotFANOVA(x = x.CW[[2]], group.label = as.character(group.label.CW))
plotFANOVA(x = x.CW[[2]], group.label = as.character(group.label.CW),
           separately = TRUE)
plotFANOVA(x = x.CW[[2]], group.label = as.character(group.label.CW),
           means = TRUE)

plotFANOVA(x = x.CW[[2]], smooth = TRUE)
plotFANOVA(x = x.CW[[2]], group.label = as.character(group.label.CW),
           smooth = TRUE)
plotFANOVA(x = x.CW[[2]], group.label = as.character(group.label.CW),
           separately = TRUE, smooth = TRUE)
plotFANOVA(x = x.CW[[2]], group.label = as.character(group.label.CW),
           means = TRUE, smooth = TRUE)}
\dontshow{
plotFANOVA(x = x.gait[[1]], int = c(0.025, 0.975))
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975))
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), separately = TRUE)
plotFANOVA(x = x.gait[[1]], group.label = as.character(group.label.gait),
           int = c(0.025, 0.975), means = TRUE)}
}
\keyword{Plot}
