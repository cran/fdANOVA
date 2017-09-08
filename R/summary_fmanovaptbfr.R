summary.fmanovaptbfr = function(object, ...){
  group.label0 = unique(object$group.label)
  l = length(group.label0); n.i = numeric(l)
  for(i in 1:l) n.i[i] = sum(object$group.label == group.label0[i])
  cat("\n", "FMANOVA - Permutation Tests based on a Basis Function Representation", "\n")
  cat("\n", "Data summary", "\n", "\n")
  if(object$basis == "own"){
    cat("Number of observations =", ncol(object$own.basis[[1]]), "\n")
    cat("Number of features =", length(object$own.basis), "\n")
    cat("Number of basis functions =", nrow(object$own.basis[[1]]), "\n")
    cat("Number of groups =", length(unique(object$group.label)), "\n")
    cat("Group labels:", group.label0, "\n")
    cat("Group sizes:", n.i, "\n")
  }else{
    cat("Number of observations =", ncol(object$data[[1]]), "\n")
    cat("Number of features =", length(object$data), "\n")
    cat("Number of time points =", nrow(object$data[[1]]), "\n")
    cat("Number of groups =", length(unique(object$group.label)), "\n")
    cat("Group labels:", group.label0, "\n")
    cat("Group sizes:", n.i, "\n")
    if(is.null(object$int)){
      cat("Range of data = [", 0, ",", nrow(object$data[[1]]), "]", "\n")
    }else{
      cat("Range of data = [", object$int[1], ",", object$int[2], "]", "\n")
    }
  }
  cat("\n", "Testing results", "\n", "\n")
  cat("W  =", object$W, " p-value =", object$pvalueW, "\n")
  cat("LH =", object$LH, " p-value =", object$pvalueLH, "\n")
  cat("P  =", object$P, " p-value =", object$pvalueP, "\n")
  cat("R  =", object$R, " p-value =", object$pvalueR, "\n")
  cat("\n", "Parameters of test", "\n", "\n")
  cat("Number of permutations =", object$B, "\n")
  if(object$basis == "own"){
    cat("Basis:", object$basis, "\n")
  }else{
    if(object$criterion == "NO"){
      if(object$basis == "Fourier"){
        cat("Basis:", object$basis, "\n")
        cat("Criterion:", object$criterion, "\n")
        cat("Method:", object$method, "\n")
        cat("Km =", object$Km, "KM =", object$KM, "maxK =", object$maxK, "\n")
      }
      if(object$basis == "b-spline"){
        cat("Basis:", object$basis, "(", "norder =", object$norder, ")", "\n")
        cat("Criterion:", object$criterion, "\n")
        cat("Method:", object$method, "\n")
        cat("Km =", object$Km, "KM =", object$KM, "maxK =", object$maxK, "\n")
      }
    }else{
      if(object$criterion == "eBIC"){
        if(object$basis == "Fourier"){
          cat("Basis:", object$basis, "\n")
          cat("Criterion:", object$criterion, "(", "gamma.eBIC =", object$gamma.eBIC, ")", "\n")
          cat("Method:", object$method, "\n")
          cat("Km =", object$Km, "KM =", object$KM, "minK =", object$minK, "maxK =", object$maxK, "\n")
        }
        if(object$basis == "b-spline"){
          cat("Basis:", object$basis, "(", "norder =", object$norder, ")", "\n")
          cat("Criterion:", object$criterion, "(", "gamma.eBIC =", object$gamma.eBIC, ")", "\n")
          cat("Method:", object$method, "\n")
          cat("Km =", object$Km, "KM =", object$KM, "minK =", object$minK, "maxK =", object$maxK, "\n")
        }
      }else{
        if(object$basis == "Fourier"){
          cat("Basis:", object$basis, "\n")
          cat("Criterion:", object$criterion, "\n")
          cat("Method:", object$method, "\n")
          cat("Km =", object$Km, "KM =", object$KM, "minK =", object$minK, "maxK =", object$maxK, "\n")
        }
        if(object$basis == "b-spline"){
          cat("Basis:", object$basis, "(", "norder =", object$norder, ")", "\n")
          cat("Criterion:", object$criterion, "\n")
          cat("Method:", object$method, "\n")
          cat("Km =", object$Km, "KM =", object$KM, "minK =", object$minK, "maxK =", object$maxK, "\n")
        }
      }
    }
  }
}
