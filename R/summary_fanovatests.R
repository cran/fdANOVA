summary.fanovatests = function(object, ...){
  group.label0 = unique(object$group.label)
  l = length(group.label0); n.i = numeric(l)
  for(i in 1:l) n.i[i] = sum(object$group.label == group.label0[i])
  cat("     Analysis of Variance for Functional Data", "\n")
  cat("\n", "Data summary", "\n", "\n")
  if(is.null(object$data)){
    cat("Number of observations =", ncol(object$FP$own.basis), "\n")
    cat("Number of basis functions =", nrow(object$FP$own.basis), "\n")
  }else{
    cat("Number of observations =", ncol(object$data), "\n")
    cat("Number of time points =", nrow(object$data), "\n")
  }
  cat("Number of groups =", length(unique(object$group.label)), "\n")
  cat("Group labels:", group.label0, "\n")
  cat("Group sizes:", n.i, "\n")
  if(!is.null(object$FP)){
    if(object$FP$basis != "own"){
      if(is.null(object$FP$int)){
        cat("Range of data = [", 0, ",", nrow(object$data), "]", "\n")
      }else{
        cat("Range of data = [", object$FP$int[1], ",", object$FP$int[2], "]", "\n")
      }
    }
  }
  cat("\n", "Testing results and parameters of tests", "\n", "\n")
  if(!is.null(object$FP)){
    cat("FP test - permutation test based on a basis function representation", "\n")
    cat("Test statistic =", object$FP$statFP, " p-value =", object$FP$pvalueFP, "\n")
    cat("Number of permutations =", object$FP$B, "\n")
    if(object$FP$basis == "own"){
      cat("Basis:", object$FP$basis, "\n", "\n")
    }else{
      if(object$FP$criterion == "NO"){
        if(object$FP$basis == "Fourier"){
          cat("Basis:", object$FP$basis, "\n")
          cat("Criterion:", object$FP$criterion, "\n")
          cat("CommonK:", object$FP$commonK, "\n")
          cat("K =", object$FP$K, "maxK =", object$FP$maxK, "\n", "\n")
        }
        if(object$FP$basis == "b-spline"){
          cat("Basis:", object$FP$basis, "(", "norder =", object$FP$norder, ")", "\n")
          cat("Criterion:", object$FP$criterion, "\n")
          cat("CommonK:", object$FP$commonK, "\n")
          cat("K =", object$FP$K, "maxK =", object$FP$maxK, "\n", "\n")
        }
      }else{
        if(object$FP$criterion == "eBIC"){
          if(object$FP$basis == "Fourier"){
            cat("Basis:", object$FP$basis, "\n")
            cat("Criterion:", object$FP$criterion, "(", "gamma.eBIC =", object$FP$gamma.eBIC, ")", "\n")
            cat("CommonK:", object$FP$commonK, "\n")
            cat("K =", object$FP$K, "minK =", object$FP$minK, "maxK =", object$FP$maxK, "\n", "\n")
          }
          if(object$FP$basis == "b-spline"){
            cat("Basis:", object$FP$basis, "(", "norder =", object$FP$norder, ")", "\n")
            cat("Criterion:", object$FP$criterion, "(", "gamma.eBIC =", object$FP$gamma.eBIC, ")", "\n")
            cat("CommonK:", object$FP$commonK, "\n")
            cat("K =", object$FP$K, "minK =", object$FP$minK, "maxK =", object$FP$maxK, "\n", "\n")
          }
        }else{
          if(object$FP$basis == "Fourier"){
            cat("Basis:", object$FP$basis, "\n")
            cat("Criterion:", object$FP$criterion, "\n")
            cat("CommonK:", object$FP$commonK, "\n")
            cat("K =", object$FP$K, "minK =", object$FP$minK, "maxK =", object$FP$maxK, "\n", "\n")
          }
          if(object$FP$basis == "b-spline"){
            cat("Basis:", object$FP$basis, "(", "norder =", object$FP$norder, ")", "\n")
            cat("Criterion:", object$FP$criterion, "\n")
            cat("CommonK:", object$FP$commonK, "\n")
            cat("K =", object$FP$K, "minK =", object$FP$minK, "maxK =", object$FP$maxK, "\n", "\n")
          }
        }
      }
    }
  }
  if(!is.null(object$CH)){
    cat("CH test - L2-norm-based parametric bootstrap test for homoscedastic samples", "\n")
    cat("Test statistic =", object$CH$statCH, " p-value =", object$CH$pvalueCH, "\n")
    cat("Number of discretized artificial trajectories for each process Z_i(t) =", object$CH$paramCH, "\n", "\n")
  }
  if(!is.null(object$CS)){
    cat("CS test - L2-norm-based parametric bootstrap test for heteroscedastic samples", "\n")
    cat("Test statistic =", object$CS$statCS, " p-value =", object$CS$pvalueCS, "\n")
    cat("Number of discretized artificial trajectories for each process Z_i(t) =", object$CS$paramCS, "\n", "\n")
  }
  if(!is.null(object$L2N)){
    cat("L2N test - L2-norm-based test with naive method of estimation", "\n")
    cat("Test statistic =", object$L2N$statL2, " p-value =", object$L2N$pvalueL2N, "\n")
    cat("beta =", object$L2N$betaL2N, " d =", object$L2N$dL2N, "\n", "\n")
  }
  if(!is.null(object$L2B)){
    cat("L2B test - L2-norm-based test with bias-reduced method of estimation", "\n")
    cat("Test statistic =", object$L2B$statL2, " p-value =", object$L2B$pvalueL2B, "\n")
    cat("beta =", object$L2B$betaL2B, " d =", object$L2B$dL2B, "\n", "\n")
  }
  if(!is.null(object$L2b)){
    cat("L2b test - L2-norm-based bootstrap test", "\n")
    cat("Test statistic =", object$L2b$statL2, " p-value =", object$L2b$pvalueL2b, "\n")
    cat("Number of bootstrap samples =", object$L2b$paramL2b, "\n", "\n")
  }
  if(!is.null(object$FN)){
    cat("FN test - F-type test with naive method of estimation", "\n")
    cat("Test statistic =", object$FN$statF, " p-value =", object$FN$pvalueFN, "\n")
    cat("d1 =", object$FN$d1FN, " d2 =", object$FN$d2FN, "\n", "\n")
  }
  if(!is.null(object$FB)){
    cat("FB test - F-type test with bias-reduced method of estimation", "\n")
    cat("Test statistic =", object$FB$statF, " p-value =", object$FB$pvalueFB, "\n")
    cat("d1 =", object$FB$d1FB, " d2 =", object$FB$d2FB, "\n", "\n")
  }
  if(!is.null(object$Fb)){
    cat("Fb test - F-type bootstrap test", "\n")
    cat("Test statistic =", object$Fb$statF, " p-value =", object$Fb$pvalueFb, "\n")
    cat("Number of bootstrap samples =", object$Fb$paramFb, "\n", "\n")
  }
  if(!is.null(object$GPF)){
    cat("GPF test - globalizing the pointwise F-test", "\n")
    cat("Test statistic =", object$GPF$statGPF, " p-value =", object$GPF$pvalueGPF, "\n")
    cat("beta =", object$GPF$betaGPF, " d =", object$GPF$dGPF, "\n", "\n")
  }
  if(!is.null(object$Fmaxb)){
    cat("Fmaxb test - Fmax bootstrap test", "\n")
    cat("Test statistic =", object$Fmaxb$statFmax, " p-value =", object$Fmaxb$pvalueFmaxb, "\n")
    cat("Number of bootstrap samples =", object$Fmaxb$paramFmaxb, "\n", "\n")
  }
  if(!is.null(object$TRP)){
    iik = 0
    for(ik in object$TRP$k){
      iik = iik + 1
      if(object$TRP$permutation == FALSE) {
        cat("TRP - tests based on k =", ik, "random projections", "\n")
        cat("p-value ANOVA =", (object$TRP$pvalues.anova)[iik], "(without permutation)", "\n")
        cat("p-value ATS =", (object$TRP$pvalues.ATS)[iik], "(without permutation)", "\n")
        cat("p-value WTPS =", (object$TRP$pvalues.WTPS)[iik], "(using B =", object$TRP$B, "permutations)", "\n", "\n")
      }else{
        cat("TRP - tests based on k =", ik, "random projections", "\n")
        cat("p-value ANOVA =", (object$TRP$pvalues.anova)[iik], "\n")
        cat("p-value ATS =", (object$TRP$pvalues.ATS)[iik], "\n")
        cat("p-value WTPS =", (object$TRP$pvalues.WTPS)[iik], "\n", "\n")
      }
    }
    if(object$TRP$permutation == TRUE){
      cat("Number of permutations =", object$TRP$B, "\n")
    }
    if(object$TRP$projection == "GAUSS"){
      cat("For projections, the Gaussian white noise was generated.", "\n")
    }else{
      cat("For projections, the Brownian motion was generated.", "\n")
    }
    if (length(object$TRP$k) > 1) {
      if (object$TRP$independent.projection.tests) {
        cat("The random projections are generated independently for different elements of vector k.", "\n")
      } else {
        cat("The random projections are generated dependently for different elements of vector k.", "\n")
      }
    }
  }
}
