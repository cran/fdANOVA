if(getRversion() >= "2.15.1")  utils::globalVariables(c("k", "test"))

plot.fanovatests = function(x, y, ...){
  if(!("ggplot2" %in% rownames(installed.packages()))){
    stop("Please install package 'ggplot2'")
  }
  if(missing(y)){
    if(class(x) != "fanovatests"){ stop("argument x is not of fanovatests class") }
    if(is.null(x$TRP)){ stop("the (standard) tests based on random projections were not performed") }
    pvalues = c(x$TRP$pvalues.anova, x$TRP$pvalues.ATS, x$TRP$pvalues.WTPS)
    da = data.frame(pvalues = pvalues, k = rep(x$TRP$k, 3),
                    test = rep(c("ANOVA", "ATS", "WTPS"), each = length(x$TRP$k)))
    figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                            geom = c("line", "point"), ylab = "p-value") +
      ggplot2::scale_shape_manual(values = c(1, 3, 5)) +
      ggplot2::scale_color_manual(values = c(1, 3, 6)) +
      ggplot2::labs(title = "FANOVA - Tests based on k Random Projections (ANOVA and ATS without permutation)")
  }else{
    if(class(y) != "fanovatests"){ stop("argument y is not of fanovatests class") }
    if(is.null(y$TRP)){ stop("the (permutation) tests based on random projections were not performed") }
    if(missing(x)){
      pvalues = c(y$TRP$pvalues.anova, y$TRP$pvalues.ATS, y$TRP$pvalues.WTPS)
      da = data.frame(pvalues = pvalues, k = rep(y$TRP$k, 3),
                      test = rep(c("ANOVAp", "ATSp", "WTPS"), each = length(y$TRP$k)))
      figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                              geom = c("line", "point"), ylab = "p-value") +
        ggplot2::scale_shape_manual(values = c(2, 4, 5)) +
        ggplot2::scale_color_manual(values = c(2, 4, 6)) +
        ggplot2::labs(title = "FANOVA - Tests based on k Random Projections (using permutations)")
    }else{
      if(class(x) != "fanovatests"){ stop("argument x is not of fanovatests class") }
      if(is.null(x$TRP)){ stop("the (standard) tests based on random projections were not performed") }
      if(any(x$TRP$k != y$TRP$k)){ stop("numbers of projections must be the same for standard and permutation tests") }
      pvalues.stand = c(x$TRP$pvalues.anova, x$TRP$pvalues.ATS, x$TRP$pvalues.WTPS)
      pvalues.perm = c(y$TRP$pvalues.anova, y$TRP$pvalues.ATS)
      da = data.frame(pvalues = c(pvalues.stand, pvalues.perm), k = rep(x$TRP$k, 5),
                      test = rep(c("ANOVA", "ATS", "WTPS", "ANOVAp", "ATSp"), each = length(x$TRP$k)))
      figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                              geom = c("line", "point"), ylab = "p-value") +
        ggplot2::scale_shape_manual(values = 1:5) +
        ggplot2::scale_color_manual(values = c(1, 2, 3, 4, 6)) +
        ggplot2::labs(title = "FANOVA - Tests based on k Random Projections")
    }
  }
  figure
}
