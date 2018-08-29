summary.fmanovatrp <- function(object, ...){
  group.label0 <- unique(object$group.label)
  l <- length(group.label0); n.i <- numeric(l)
  for(i in 1:l) n.i[i] <- sum(object$group.label == group.label0[i])
  cat("\n", "FMANOVA - Tests based on Random Projections", "\n")
  cat("\n", "Data summary", "\n", "\n")
  cat("Number of observations =", ncol(object$data[[1]]), "\n")
  cat("Number of features =", length(object$data), "\n")
  cat("Number of time points =", nrow(object$data[[1]]), "\n")
  cat("Number of groups =", length(unique(object$group.label)), "\n")
  cat("Group labels:", group.label0, "\n")
  cat("Group sizes:", n.i, "\n")
  cat("\n", "Testing results", "\n")
  iik <- 0
  for(ik in object$k){
    iik <- iik + 1
    cat("\n", "Number of random projections k =", ik, "\n")
    cat("p-value Wilks =", (object$pvalues)[1, iik], "\n")
    cat("p-value Lawley-Hotelling =", (object$pvalues)[2, iik], "\n")
    cat("p-value Pillai =", (object$pvalues)[3, iik], "\n")
    cat("p-value Roy =", (object$pvalues)[4, iik], "\n")
  }
  cat("\n", "Parameters of test", "\n", "\n")
  if(object$permutation == TRUE){
    cat("Number of permutations =", object$B, "\n")
  }
  if(object$projection == "GAUSS"){
    cat("For projections, the Gaussian white noise was generated.", "\n")
  }else{
    cat("For projections, the Brownian motion was generated.", "\n")
  }
  if (length(object$k) > 1) {
    if (object$independent.projection.tests) {
      cat("The random projections are generated independently for different elements of vector k.", "\n")
    } else {
      cat("The random projections are generated dependently for different elements of vector k.", "\n")
    }
  }
}
