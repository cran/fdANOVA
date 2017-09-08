print.fmanovaptbfr = function(x, ...){
  cat("\n", "FMANOVA - Permutation Tests based on a Basis Function Representation", "\n", "\n")
  cat("W  =", x$W, " p-value =", x$pvalueW, "\n")
  cat("LH =", x$LH, " p-value =", x$pvalueLH, "\n")
  cat("P  =", x$P, " p-value =", x$pvalueP, "\n")
  cat("R  =", x$R, " p-value =", x$pvalueR, "\n")
}
