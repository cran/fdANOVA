if(getRversion() >= "2.15.1")  utils::globalVariables(c("values1", "a", "values", "t1"))

plotFANOVA = function(x, group.label = NULL, int = NULL, separately = FALSE, means = FALSE, smooth = FALSE, ...){
  if(!("ggplot2" %in% rownames(installed.packages()))){
    stop("Please install package 'ggplot2'")
  }
  group.label0 = unique(group.label)
  l = length(group.label0)
  if(is.null(int)){
    tt = seq_len(nrow(x))
  }else{
    if(length(int) != 2){ stop("argument int must be of length two") }
    tt = seq(int[1], int[2], length = nrow(x))
  }
  if(smooth){
    x.ss = matrix(0, nrow = nrow(x), ncol = ncol(x))
    for(ii in 1:ncol(x)){
      ss = smooth.spline(tt, x[, ii])
      x.ss[, ii] = predict(ss, tt)$y
    }
    x = x.ss
  }
  if(means){
    d1 = data.frame(t = rep(tt, ncol(x)))
    d1$a = rep(seq_len(ncol(x)), each = nrow(x))
    d1$group.label = rep(group.label, each = nrow(x))
    d1$values1 = rep(0, nrow(d1))
    for(i in 1:l){
      d1$values1[d1$group.label == group.label0[i]] = rowMeans(x[, group.label == group.label0[i]])
    }
    ggplot2::qplot(t, values1, data = d1, group = a, geom = "line", ylab = '', color = group.label, size = I(1.5))
  }else{
    data2 = stack(as.data.frame(x))
    data2$a = rep(seq_len(ncol(x)), each = nrow(x))
    data2$group.label = rep(group.label, each = nrow(x))
    data2$t = rep(tt, ncol(x))
    if(separately){
      if(is.null(group.label)){ stop("invalid argument group.label") }
      d1 = data.frame(t1 = data2$t)
      d1$a = rep(seq_len(ncol(x)), each = nrow(x))
      d1$group.label = rep(group.label, each = nrow(x))
      d1$values1 = rep(0, nrow(d1))
      for(i in 1:l){
        d1$values1[d1$group.label == group.label0[i]] = rowMeans(x[, group.label == group.label0[i]])
      }
      ggplot2::qplot(t, values, data = data2, group = a, geom = "line", ylab = '', facets = group.label ~ ., colour = I("cornflowerblue")) +
        ggplot2::geom_line(ggplot2::aes(t1, values1, group = a), data = d1, colour = "tomato", size = I(1.5))
    }else{
      if(is.null(group.label)){
        ggplot2::qplot(t, values, data = data2, group = a, geom = "line", ylab = '') +
          ggplot2::geom_line(size = I(0.000005), colour = I("cornflowerblue"))
      }else{
        d1 = data.frame(t1 = data2$t)
        d1$a = rep(seq_len(ncol(x)), each = nrow(x))
        d1$group.label = rep(group.label, each = nrow(x))
        d1$values1 = rep(0, nrow(d1))
        for(i in 1:l){
          d1$values1[d1$group.label == group.label0[i]] = rowMeans(x[, group.label == group.label0[i]])
        }
        ggplot2::qplot(t, values, data = data2, group = a, geom = "line", ylab = '', color = group.label) +
          ggplot2::geom_line(ggplot2::aes(t1, values1, group = a), data = d1, size = I(1.5))
      }
    }
  }
}
