print.fanovatests = function(x, ...){
  cat("     Analysis of Variance for Functional Data", "\n", "\n")
  if(!is.null(x$FP)){
    cat("FP test - permutation test based on a basis function representation", "\n")
    cat("Test statistic =", x$FP$statFP, " p-value =", x$FP$pvalueFP, "\n", "\n")
  }
  if(!is.null(x$CH)){
    cat("CH test - L2-norm-based parametric bootstrap test for homoscedastic samples", "\n")
    cat("Test statistic =", x$CH$statCH, " p-value =", x$CH$pvalueCH, "\n", "\n")
  }
  if(!is.null(x$CS)){
    cat("CS test - L2-norm-based parametric bootstrap test for heteroscedastic samples", "\n")
    cat("Test statistic =", x$CS$statCS, " p-value =", x$CS$pvalueCS, "\n", "\n")
  }
  if(!is.null(x$L2N)){
    cat("L2N test - L2-norm-based test with naive method of estimation", "\n")
    cat("Test statistic =", x$L2N$statL2, " p-value =", x$L2N$pvalueL2N, "\n", "\n")
  }
  if(!is.null(x$L2B)){
    cat("L2B test - L2-norm-based test with bias-reduced method of estimation", "\n")
    cat("Test statistic =", x$L2B$statL2, " p-value =", x$L2B$pvalueL2B, "\n", "\n")
  }
  if(!is.null(x$L2b)){
    cat("L2b test - L2-norm-based bootstrap test", "\n")
    cat("Test statistic =", x$L2b$statL2, " p-value =", x$L2b$pvalueL2b, "\n", "\n")
  }
  if(!is.null(x$FN)){
    cat("FN test - F-type test with naive method of estimation", "\n")
    cat("Test statistic =", x$FN$statF, " p-value =", x$FN$pvalueFN, "\n", "\n")
  }
  if(!is.null(x$FB)){
    cat("FB test - F-type test with bias-reduced method of estimation", "\n")
    cat("Test statistic =", x$FB$statF, " p-value =", x$FB$pvalueFB, "\n", "\n")
  }
  if(!is.null(x$Fb)){
    cat("Fb test - F-type bootstrap test", "\n")
    cat("Test statistic =", x$Fb$statF, " p-value =", x$Fb$pvalueFb, "\n", "\n")
  }
  if(!is.null(x$GPF)){
    cat("GPF test - globalizing the pointwise F-test", "\n")
    cat("Test statistic =", x$GPF$statGPF, " p-value =", x$GPF$pvalueGPF, "\n", "\n")
  }
  if(!is.null(x$Fmaxb)){
    cat("Fmaxb test - Fmax bootstrap test", "\n")
    cat("Test statistic =", x$Fmaxb$statFmax, " p-value =", x$Fmaxb$pvalueFmaxb, "\n", "\n")
  }
  if(!is.null(x$TRP)){
    iik = 0
    for(ik in x$TRP$k){
      iik = iik + 1
      cat("TRP - tests based on k =", ik, "random projections", "\n")
      cat("p-value ANOVA =", (x$TRP$pvalues.anova)[iik], "\n")
      cat("p-value ATS =", (x$TRP$pvalues.ATS)[iik], "\n")
      cat("p-value WTPS =", (x$TRP$pvalues.WTPS)[iik], "\n", "\n")
    }
  }
}
