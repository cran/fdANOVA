print.fmanovatrp = function(x, ...){
  cat("\n", "FMANOVA - Tests based on Random Projections", "\n")
  iik = 0
  for(ik in x$k){
    iik = iik + 1
    cat("\n", "Number of random projections k =", ik, "\n")
    cat("p-value Wilks =", (x$pvalues)[1, iik], "\n")
    cat("p-value Lawley-Hotelling =", (x$pvalues)[2, iik], "\n")
    cat("p-value Pillai =", (x$pvalues)[3, iik], "\n")
    cat("p-value Roy =", (x$pvalues)[4, iik], "\n")
  }
}
