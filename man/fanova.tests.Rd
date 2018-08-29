\name{fanova.tests}
\alias{fanova.tests}
\title{
Tests for FANOVA Problem
}
\description{
Performs the testing procedures for the one-way analysis of variance for (univariate) functional data (FANOVA). See Section 2.1 of the vignette file (\code{vignette("fdANOVA", package = "fdANOVA")}), for details of the tests.

We consider the \eqn{l} groups of independent random functions \eqn{X_{ij}(t)}, \eqn{i=1,\dots,l,} \eqn{j=1,\dots,n_i} defined over a closed and bounded interval \eqn{I=[a,b]}. Let \eqn{n=n_1+\dots+n_l}. These groups may differ in mean functions, i.e., we assume that \eqn{X_{ij}(t)}, \eqn{j=1,\dots,n_i} are stochastic processes with mean function \eqn{\mu_i(t)}, \eqn{t\in I} and covariance function \eqn{\gamma(s, t)}, \eqn{s,t\in I}, for \eqn{i=1,\dots,l}. Of interest is to test the following null hypothesis
\deqn{
H_0:\mu_1(t)=\dots=\mu_l(t),\ t\in I.
}
The alternative is the negation of the null hypothesis. We assume that each functional observation is observed on a common grid of \eqn{\mathcal{T}} design time points equally spaced in \eqn{I} (see Section 3.1 of the vignette file, \code{vignette("fdANOVA", package = "fdANOVA")}).
}
\usage{
fanova.tests(x = NULL, group.label, test = "ALL",
             params = NULL,
             parallel = FALSE, nslaves = NULL)

# more detailed usage of params:
# params = list(paramFP = list(int, B.FP = 1000,
#                              basis = c("Fourier", "b-spline", "own"),
#                              own.basis, own.cross.prod.mat,
#                              criterion = c("BIC", "eBIC", "AIC", "AICc", "NO"),
#                              commonK = c("mode", "min", "max", "mean"),
#                              minK = NULL, maxK = NULL, norder = 4, gamma.eBIC = 0.5)
#               paramCH = 10000,
#               paramCS = 10000,
#               paramL2b = 10000,
#               paramFb = 10000,
#               paramFmaxb = 10000,
#               paramTRP = list(k = 30, projection = c("GAUSS", "BM"),
#                               permutation = FALSE, B.TRP = 10000,
#                               independent.projection.tests = TRUE))
}
\arguments{
\item{x}{a \eqn{\mathcal{T}\times n} matrix of data, whose each column is a discretized version of a function and rows correspond to design time points. Its default values is \code{NULL}, since if the FP test is only used, we can give a basis representation of the data instead of raw observations (see the list \code{paramFP} below). For any of the other testing procedures, the raw data are needed.
}
\item{group.label}{a vector containing group labels.
}
\item{test}{a kind of indicator which establishes a choice of FANOVA tests to be performed. Its default value means that all testing procedures of Section 2.1 of the vignette file will be used. When we want to use only some tests, the parameter \code{test} is an appropriate subvector of the following vector of tests' labels \code{c("FP",} \code{"CH",} \code{"CS",} \code{"L2N",} \code{"L2B",} \code{"L2b",} \code{"FN",} \code{"FB",} \code{"Fb",} \code{"GPF",} \code{"Fmaxb",} \code{"TRP")}, where \code{"FP"} - permutation test based on basis function representation (Gorecki and Smaga, 2015); \code{"CH"} and \code{"CS"} - L2-norm-based parametric bootstrap tests for homoscedastic and heteroscedastic samples, respectively (Cuevas et al., 2004); \code{"L2N"} and \code{"L2B"} - L2-norm-based test with naive and bias-reduced method of estimation, respectively (Faraway, 1997; Zhang and Chen, 2007; Zhang, 2013); \code{"L2b"} - L2-norm-based bootstrap test (Zhang, 2013); \code{"FN"} and \code{"FB"} - F-type test with naive and bias-reduced method of estimation, respectively (Shen and Faraway, 2004; Zhang, 2011); \code{"Fb"} - F-type bootstrap test (Zhang, 2013); \code{"GPF"} - globalizing the pointwise F-test (Zhang and Liang, 2014); \code{"Fmaxb"} - Fmax bootstrap test (Zhang et al., 2018); \code{"TRP"} - tests based on random projections (Cuesta-Albertos and Febrero-Bande, 2010).
}
\item{params}{a list of additional parameters for the FP, CH, CS, L\eqn{^2}b, Fb, Fmaxb tests and the tests based on random projections. It can contain all or a part of the elements \code{paramFP}, \code{paramCH}, \code{paramCS}, \code{paramL2b}, \code{paramFb}, \code{paramFmaxb} and \code{paramTRP} for passing the parameters for the FP, CH, CS, L\eqn{^2}b, Fb, Fmaxb tests and tests based on random projections, respectively, to the function \code{fanova.tests}. They are described below. The default value of \code{params} means that these tests are performed with their default values.
}
\item{paramFP}{a list containing the parameters for the FP test.}
\item{int}{a vector of two elements representing the interval \eqn{I=[a,b]}. When it is not specified, it is determined by a number of design time points.
}
\item{B.FP}{a number of permutation replicates for the FP tests.}
\item{basis}{a choice of basis of functions used in the basis function representation of the data.
}
\item{own.basis}{if \code{basis = "own"}, a \eqn{K\times n} matrix with columns containing the coefficients of the basis function representation of the observations.
}
\item{own.cross.prod.mat}{if \code{basis = "own"}, a \eqn{K\times K} cross product matrix corresponding to a basis used to obtain the matrix \code{own.basis}.
}
\item{criterion}{a choice of information criterion for selecting the optimum value of \eqn{K}.

\code{criterion = "NO"} means that \eqn{K} is equal to the parameter \code{maxK} defined below. We have \deqn{\code{BIC}(X_{ij})=\mathcal{T}\log(\mathbf{e}_{ij}^{\top}\mathbf{e}_{ij}/\mathcal{T})+K\log\mathcal{T},} \deqn{\code{eBIC}(X_{ij})=\mathcal{T}\log(\mathbf{e}_{ij}^{\top}\mathbf{e}_{ij}/\mathcal{T})+K[\log\mathcal{T}+2\gamma\log(K_{\max})],} \deqn{\code{AIC}(X_{ij})=\mathcal{T}\log(\mathbf{e}_{ij}^{\top}\mathbf{e}_{ij}/\mathcal{T})+2K} and \deqn{\code{AICc}(X_{ij})=\code{AIC}(X_{ij})+2K(K + 1)/(n-K-1),} where \deqn{\mathbf{e}_{ij}=(e_{ij1},\dots,e_{ij\mathcal{T}})^{\top},} \deqn{e_{ijr}=X_{ij}(t_r)-\sum_{m=1}^K\hat{c}_{ijm}\varphi_m(t_r),} \eqn{t_1,\dots,t_{\mathcal{T}}} are the design time points, \eqn{\gamma\in[0,1]}, \eqn{K_{\max}} is a maximum \eqn{K} considered and \eqn{\log} denotes the natural logarithm.
}
\item{commonK}{a choice of method for selecting the common value for all observations from the values of \eqn{K} corresponding to all processes.
}
\item{minK}{a minimum value of \eqn{K}. When \code{basis = "Fourier"}, it has to be an odd number. If \code{minK = NULL}, we take \code{minK = 3}. For \code{basis = "b-spline"}, \code{minK} has to be greater than or equal to \code{norder} defined below. If \code{minK = NULL} or \code{minK < norder}, then we take \code{minK = norder}.
}
\item{maxK}{a maximum value of \eqn{K}. When \code{basis = "Fourier"}, it has to be an odd number. If \code{maxK = NULL}, we take \code{maxK} equal to the largest odd number smaller than the number of design time points. If \code{maxK} is greater than or equal to the number of design time points, \code{maxK} is taken as above. For \code{basis = "b-spline"}, \code{maxK} has to be smaller than or equal to the number of design time points. If \code{maxK = NULL} or \code{maxK} is greater than the number of design time points, then we take \code{maxK} equal to the number of design time points.
}
\item{norder}{if \code{basis = "b-spline"}, an integer specifying the order of b-splines.}
\item{gamma.eBIC}{a \eqn{\gamma\in[0,1]} parameter in the eBIC.}
\item{paramCH}{a number of discretized artificial trajectories for generating Gaussian processes for the CH test.}
\item{paramCS}{a number of discretized artificial trajectories for generating Gaussian processes for the CS test.}
\item{paramL2b}{a number of bootstrap samples for the L\eqn{^2}b test.}
\item{paramFb}{a number of bootstrap samples for the Fb test.}
\item{paramFmaxb}{a number of bootstrap samples for the Fmaxb test.}
\item{paramTRP}{a list containing the parameters of the tests based on random projections.}
\item{k}{a vector of numbers of projections.}
\item{projection}{a method of generating Gaussian processes in step 1 of the tests based on random projections presented in Section 2 of the vignette file. If \code{projection = "GAUSS"}, the Gaussian white noise is generated as in the function \code{anova.RPm} from the R package \pkg{fda.usc}. In the second case, the Brownian motion is generated.}
\item{permutation}{a logical indicating whether to compute p-values of the tests based on random projections by permutation method.}
\item{B.TRP}{a number of permutation replicates for the tests based on random projections.}
\item{independent.projection.tests}{a logical indicating whether to generate the random projections independently or dependently for different elements of vector \code{k}. In the first case, the random projections for each element of vector \code{k} are generated separately, while in the second one, they are generated as chained subsets, e.g., for \code{k = c(5, 10)}, the first 5 projections are a subset of the last 10. The second way of generating random projections is faster than the first one.}
\item{parallel}{a logical indicating whether to use parallelization.}
\item{nslaves}{if \code{parallel = TRUE}, a number of slaves. Its default value means that it will be equal to a number of logical processes of a computer used.
}
}
\details{
To perform step 3 of the projection procedure given in Section 2.1 of the vignette file, we use five tests: the standard (\code{paramTRP$permutation = FALSE}) and permutation (\code{paramTRP$permutation} \code{= TRUE}) tests based on ANOVA F-test statistic and ANOVA-type statistic (ATS) proposed by Brunner et al. (1997), as well as the testing procedure based on Wald-type permutation statistic (WTPS) of Pauly et al. (2015).
}
\value{
A list with class "fanovatests" containing the following components (|\code{k}| denotes the length of vector \code{k}):
\item{FP}{a list containing value of test statistic \code{statFP}, p-value \code{pvalueFP} and used parameters for the FP test. The chosen optimal length of basis expansion \code{K} is also given there.
}
\item{CH}{a list containing value of test statistic \code{statCH}, p-value \code{pvalueCH} and used parameter \code{paramCH} for the CH test.
}
\item{CS}{a list containing value of test statistic \code{statCS}, p-value \code{pvalueCS} and used parameter \code{paramCS} for the CS test.
}
\item{L2N}{a list containing value of test statistic \code{statL2}, p-value \code{pvalueL2N} and values of estimators \code{betaL2N} and \code{dL2N} used in approximation of null distribution of test statistic for the L\eqn{^2}N test.
}
\item{L2B}{a list containing value of test statistic \code{statL2}, p-value \code{pvalueL2B} and values of estimators \code{betaL2B} and \code{dL2B} used in approximation of null distribution of test statistic for the L\eqn{^2}B test.
}
\item{L2b}{a list containing value of test statistic \code{statL2}, p-value \code{pvalueL2b} and used parameter \code{paramL2b} for the L\eqn{^2}b test.
}
\item{FN}{a list containing value of test statistic \code{statF}, p-value \code{pvalueFN} and values of estimators \code{d1FN} and \code{d2FN} used in approximation of null distribution of test statistic for the FN test.
}
\item{FB}{a list containing value of test statistic \code{statF}, p-value \code{pvalueFB} and values of estimators \code{d1FB} and \code{d2FB} used in approximation of null distribution of test statistic for the FB test.
}
\item{Fb}{a list containing value of test statistic \code{statF}, p-value \code{pvalueFb} and used parameter \code{paramFb} for the Fb test.
}
\item{GPF}{a list containing value of test statistic \code{statGPF}, p-value \code{pvalueGPF} and values of estimators \code{betaGPF} and \code{dGPF} used in approximation of null distribution of test statistic for the GPF test.
}
\item{Fmaxb}{a list containing value of test statistic \code{statFmax}, p-value \code{pvalueFmaxb} and used parameter \code{paramFmaxb} for the Fmaxb test.
}
\item{TRP}{a list containing the following elements: vectors \code{pvalues.anova}, \code{pvalues.ATS}, \code{pvalues.WTPS} of length |\code{k}| containing p-values for tests based on random projections and for numbers of projections given in \code{k}; if \code{independent.projection.tests = TRUE}, a list \code{data.projections} of length |\code{k}|, whose \eqn{i}th element is an \eqn{n\times} \code{k[i]} matrix with columns being projections of the data; when \code{independent.projection.tests = FALSE}, an \eqn{n\times \max}(\code{k}) matrix \code{data.projections} with columns being projections of the data; used parameters for the tests based on random projections.
}
and the values of other used parameters: \code{data = x}, \code{group.label}, etc.
}
\references{
Brunner E, Dette H, Munk A (1997). Box-Type Approximations in Nonparametric Factorial
Designs. Journal of the American Statistical Association 92, 1494-1502.

Cuesta-Albertos JA, Febrero-Bande M (2010). A Simple Multiway ANOVA for Functional
Data. Test 19, 537-557.

Cuevas A, Febrero M, Fraiman R (2004). An Anova Test for Functional Data. Computational
Statistics & Data Analysis 47, 111-122.

Faraway J (1997). Regression Analysis for a Functional Response. Technometrics 39,
254-261.

Gorecki T, Smaga L (2015). A Comparison of Tests for the One-Way ANOVA Problem for
Functional Data. Computational Statistics 30, 987-1010.

Gorecki T, Smaga L (2017). Multivariate Analysis of Variance for Functional Data. Journal of Applied Statistics 44, 2172-2189.

Pauly M, Brunner E, Konietschke F (2015). Asymptotic Permutation Tests in General
Factorial Designs. Journal of the Royal Statistical Society Series B 77, 461-473.

Shen Q, Faraway J (2004). An F Test for Linear Models with Functional Responses. Statistica Sinica 14, 1239-1257.

Zhang JT (2011). Statistical Inferences for Linear Models with Functional Responses. Statistica Sinica 21, 1431-1451.

Zhang JT (2013). Analysis of Variance for Functional Data. Chapman & Hall, London.

Zhang JT, Chen JW (2007). Statistical Inferences for Functional Data. The Annals of
Statistics 35, 1052-1079.

Zhang JT, Cheng MY, Wu HT, Zhou B (2018). A New Test for Functional One-way
ANOVA with Applications to Ischemic Heart Screening. Computational Statistics and Data Analysis
https://doi.org/10.1016/j.csda.2018.05.004

Zhang JT, Liang X (2014). One-Way ANOVA for Functional Data via Globalizing the
Pointwise F-Test. Scandinavian Journal of Statistics 41, 51-71.
}
\author{
Tomasz Gorecki, Lukasz Smaga
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fmanova.ptbfr}}, \code{\link{fmanova.trp}}, \code{\link{plotFANOVA}}, \code{\link{plot.fanovatests}}
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
# all FANOVA tests with default parameters
set.seed(123)
(fanova1 <- fanova.tests(x = x.gait, group.label = group.label.gait))
summary(fanova1)
# data projections generated in the test based on random projections
fanova1$TRP$data.projections

# only three tests with non-default parameters
set.seed(123)
fanova2 <- fanova.tests(x.gait, group.label.gait,
                        test = c("FP", "GPF", "Fmaxb"),
                        params = list(paramFP = list(int = c(0.025, 0.975),
                                                     B.FP = 1000, basis = "b-spline",
                                                     criterion = "eBIC",
                                                     commonK = "mean",
                                                     minK = 5, maxK = 20,
                                                     norder = 4, gamma.eBIC = 0.7),
                                      paramFmaxb = 1000))
summary(fanova2)

# the FP test with predefined basis function representation
library(fda)
fbasis <- create.bspline.basis(rangeval = c(0.025, 0.975), 19, norder = 4)
own.basis <- Data2fd(seq(0.025, 0.975, length = 20), x.gait, fbasis)$coefs
own.cross.prod.mat <- inprod(fbasis, fbasis)
set.seed(123)
fanova3 <- fanova.tests(group.label = group.label.gait, test = "FP",
                        params = list(paramFP = list(B.FP = 1000, basis = "own",
                                                     own.basis = own.basis,
                                                     own.cross.prod.mat = own.cross.prod.mat)))
summary(fanova3)

# the tests based on random projections with the Gaussian white noise generated for projections
set.seed(123)
fanova4 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30), B.TRP = 1000)))
summary(fanova4)
set.seed(123)
fanova5 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30),
                                                      permutation = TRUE, B.TRP = 1000)))
summary(fanova5)

# the tests based on random projections with the Brownian motion generated for projections
set.seed(123)
fanova6 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30), projection = "BM",
                                                      B.TRP = 1000)))
summary(fanova6)
set.seed(123)
fanova7 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        parallel = TRUE, nslaves = 2,
                        params = list(paramTRP = list(k = c(10, 20, 30), projection = "BM",
                                                      permutation = TRUE, B.TRP = 1000)))
summary(fanova7)}
\dontshow{
# only three tests with non-default parameters
fanova2 <- fanova.tests(x.gait, group.label.gait,
                        test = c("FP", "GPF", "Fmaxb"),
                        params = list(paramFP = list(int = c(0.025, 0.975),
                                                     B.FP = 3, basis = "b-spline",
                                                     criterion = "eBIC",
                                                     commonK = "mean",
                                                     minK = 5, maxK = 10,
                                                     norder = 4, gamma.eBIC = 0.7),
                                      paramFmaxb = 3))

# the FP test with predefined basis function representation
fbasis <- create.bspline.basis(rangeval = c(0.025, 0.975), 19, norder = 4)
own.basis <- Data2fd(seq(0.025, 0.975, length = 20), x.gait, fbasis)$coefs
own.cross.prod.mat <- inprod(fbasis, fbasis)
fanova3 <- fanova.tests(group.label = group.label.gait, test = "FP",
                        params = list(paramFP = list(B.FP = 3, basis = "own",
                                                     own.basis = own.basis,
                                                     own.cross.prod.mat = own.cross.prod.mat)))

# the tests based on random projections with the Gaussian white noise generated for projections
fanova4 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        params = list(paramTRP = list(k = c(1, 2), B.TRP = 2)))

# data projections generated in the test based on random projections
fanova4$TRP$data.projections

# the tests based on random projections with the Brownian motion generated for projections
fanova6 <- fanova.tests(x.gait, group.label.gait, test = "TRP",
                        params = list(paramTRP = list(k = c(1, 2), projection = "BM",
                                                      B.TRP = 2)))}
}
\keyword{Multivariate}
\keyword{Test}
