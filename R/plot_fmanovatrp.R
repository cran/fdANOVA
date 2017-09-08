if(getRversion() >= "2.15.1")  utils::globalVariables(c("k", "test"))

plot.fmanovatrp = function(x, y, withoutRoy = FALSE, ...){
  if(!("ggplot2" %in% rownames(installed.packages()))){
    stop("Please install package 'ggplot2'")
  }
  if(! is.logical(withoutRoy)){ stop("argument withoutRoy is not logical") }
  if(missing(y)){
    if(class(x) != "fmanovatrp"){ stop("argument x is not of fmanovatrp class") }
    if(withoutRoy == FALSE){
      pvalues = c()
      for(ii in 1:4){ pvalues = c(pvalues, (x$pvalues)[ii, ]) }
      da = data.frame(pvalues = pvalues, k = rep(x$k, 4),
                      test = rep(c("Wp", "LHp", "Pp", "Rp"), each = length(x$k)))
      figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                              geom = c("line", "point"), ylab = "p-values") +
        ggplot2::scale_shape_manual(values = c(1, 3, 5, 7)) +
        ggplot2::scale_color_manual(values = c(1, 3, 7, 6)) +
        ggplot2::labs(title = "FMANOVA - Tests based on k Random Projections (without permutation)")
    }else{
      pvalues = c()
      for(ii in 1:3){ pvalues = c(pvalues, (x$pvalues)[ii, ]) }
      da = data.frame(pvalues = pvalues, k = rep(x$k, 3),
                      test = rep(c("Wp", "LHp", "Pp"), each = length(x$k)))
      figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                              geom = c("line", "point"), ylab = "p-values") +
        ggplot2::scale_shape_manual(values = c(1, 3, 7)) +
        ggplot2::scale_color_manual(values = c(1, 3, 6)) +
        ggplot2::labs(title = "FMANOVA - Tests based on k Random Projections (without permutation)")
    }
  }else{
    if(class(y) != "fmanovatrp"){ stop("argument y is not of fmanovatrp class") }
    if(missing(x)){
      pvalues = c()
      for(ii in 1:4){ pvalues = c(pvalues, (y$pvalues)[ii, ]) }
      da = data.frame(pvalues = pvalues, k = rep(y$k, 4),
                      test = rep(c("Wpp", "LHpp", "Ppp", "Rpp"), each = length(y$k)))
      figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                              geom = c("line", "point"), ylab = "p-values") +
        ggplot2::scale_shape_manual(values = c(2, 4, 6, 8)) +
        ggplot2::scale_color_manual(values = c(2, 4, 5, 8)) +
        ggplot2::labs(title = "FMANOVA - Tests based on k Random Projections (using permutations)")
    }else{
      if(class(x) != "fmanovatrp"){ stop("argument x is not of fmanovatrp class") }
      if(any(x$k != y$k)){ stop("numbers of projections must be the same for standard and permutation tests") }
      if(withoutRoy == FALSE){
        pvalues.stand = c(); pvalues.perm = c()
        for(ii in 1:4){ pvalues.stand = c(pvalues.stand, (x$pvalues)[ii, ]); pvalues.perm = c(pvalues.perm, (y$pvalues)[ii, ]) }
        da = data.frame(pvalues = c(pvalues.stand, pvalues.perm), k = rep(x$k, 8),
                        test = rep(c("Wp", "LHp", "Pp", "Rp", "Wpp", "LHpp", "Ppp", "Rpp"), each = length(x$k)))
        figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                                geom = c("line", "point"), ylab = "p-values") +
          ggplot2::scale_shape_manual(values = 1:8) +
          ggplot2::scale_color_manual(values = c(1:4, 7, 5, 6, 8)) +
          ggplot2::labs(title = "FMANOVA - Tests based on k Random Projections")
      }else{
        pvalues.stand = c(); pvalues.perm = c()
        for(ii in 1:3){ pvalues.stand = c(pvalues.stand, (x$pvalues)[ii, ]) }
        for(ii in 1:4){ pvalues.perm = c(pvalues.perm, (y$pvalues)[ii, ]) }
        da = data.frame(pvalues = c(pvalues.stand, pvalues.perm), k = rep(x$k, 7),
                        test = rep(c("Wp", "LHp", "Pp", "Wpp", "LHpp", "Ppp", "Rpp"), each = length(x$k)))
        figure = ggplot2::qplot(k, pvalues, data = da, group = test, colour = test, shape = test,
                                geom = c("line", "point"), ylab = "p-values") +
          ggplot2::scale_shape_manual(values = c(1:4, 6:8)) +
          ggplot2::scale_color_manual(values = c(1:4, 5, 6, 8)) +
          ggplot2::labs(title = "FMANOVA - Tests based on k Random Projections")
      }
    }
  }
  figure
}
