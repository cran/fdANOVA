fanova.tests = function(x = NULL, group.label, test = "ALL", params = NULL,
                        parallel = FALSE, nslaves = NULL){
  if(any(is.na(group.label))){ stop("argument group.label can not contain NA values") }
  if(test[1] == "ALL"){
    test = c("FP", "CH", "CS", "L2N", "L2B", "L2b", "FN", "FB", "Fb", "GPF", "Fmaxb", "TRP")
  }else{
    for(i in 1:length(test)){
      if(!(test[i] %in% c("FP", "CH", "CS", "L2N", "L2B", "L2b", "FN", "FB", "Fb", "GPF", "Fmaxb", "TRP"))){
        stop("argument test must be the subvector of c('FP', 'CH', 'CS', 'L2N', 'L2B', 'L2b', 'FN', 'FB', 'Fb', 'GPF', 'Fmaxb', 'TRP')")
      }
    }
  }
  group.label0 = unique(group.label)
  l = length(group.label0)
  n.i = numeric(l)
  for(i in 1:l) n.i[i] = sum(group.label == group.label0[i])
  anova.statistic.quick = function(x, group.label){
    n = length(c(x))
    a = length(unique(group.label))
    y = data.frame(cbind(x, group.label))
    colnames(y) = c("data", "group.label")
    means = (doBy::summaryBy(data ~ group.label, data = y))[,2]
    mean.g = mean(x)
    SSA = sum(table(group.label) * ((means - mean.g)^2))
    SST = sum((x - mean.g)^2)
    SSE = SST - SSA
    return((n-a)*SSA/SSE/(a-1))
  }

  ATS.simple = function(x, group.label){
    N = length(x); n = table(group.label); ng = length(as.vector(n))
    means = matrix(tapply(x, group.label, mean), ncol = 1)
    MM = diag(1, ng) - matrix(1, ncol = ng, nrow = ng)/ng
    DM = diag(diag(MM))
    SN = N * diag(as.vector((tapply(x, group.label, var) * (n - 1)/n^2)))
    Delta = diag(as.vector(1/(n - 1)))
    test.stat = N * t(means) %*% MM %*% means/sum(diag(MM %*% SN))
    f1 = sum(diag(MM %*% SN))^2/sum(diag(MM %*% SN %*% MM %*% SN))
    f0 = sum(diag(MM %*% SN))^2/sum(diag(DM %*% DM %*% SN %*% SN %*% Delta))
    return(c(test.stat, 1 - pf(test.stat, f1, f0)))
  }

  WTPSp = function(x, group.label, group.label0, n, n.i, l, perm.WTPS, nrTRP){
    means = matrix(0, nrow = l, ncol = 1); variances = numeric(l)
    for(i in 1:l){
      means[i,] = mean(x[group.label == group.label0[i]])
      variances[i] = var(x[group.label == group.label0[i]])
    }
    V.N = n*diag(variances/n.i)
    TT = diag(rep(1, l)) - (1/l)*matrix(1, ncol = l, nrow = l)
    Q.N.T = n * t(means) %*% TT %*% MASS::ginv(TT %*% V.N %*% TT) %*% TT %*% means
    A = t(rep(1/n.i[1], n.i[1]))
    A1 = t(rep(1, n.i[1]))
    for(i in 2:l){
      A = magic::adiag(A, t(rep(1/n.i[i], n.i[i])))
      A1 = magic::adiag(A1, t(rep(1, n.i[i])))
    }
    xperm = matrix(x[perm.WTPS], ncol = nrTRP)
    meansP = A %*% xperm
    varsP = (A1 %*% (xperm^2) - n.i * meansP^2)/(n.i * (n.i - 1))
    statistics.p = sapply(1:nrTRP, function(arg){
      TP = TT %*% MASS::ginv(TT %*% (n * diag(varsP[, arg])) %*% TT) %*% TT
      statistics.p = diag(n * t(meansP[, arg]) %*% TP %*% meansP[, arg])
    })
    return(mean(statistics.p >= c(Q.N.T)))
  }
  if(any(c("FP", "CH", "CS", "L2b", "Fb", "Fmaxb", "TRP") %in% test)){
    if(!parallel){
      parallel.method = "parallel.method0"
    }else{
      if(!("doParallel" %in% rownames(installed.packages()))){
        stop("Please install package 'doParallel'")
      }
      # require(foreach, quietly = TRUE)
      requireNamespace("foreach", quietly = TRUE)
      nlp = parallel::detectCores()
      if(is.null(nslaves)){
        if(nlp >= 2){
          nslaves = nlp
          parallel.method = "parallel.method1"
        }else{
          parallel.method = "parallel.method0"
        }
      }else{
        if(nlp >= 2){
          if(nslaves >= 2){
            nslaves = nslaves
          }else{
            nslaves = nlp
          }
          parallel.method = "parallel.method1"
        }else{
          parallel.method = "parallel.method0"
        }
      }
      if(parallel.method == "parallel.method1"){
        cl = parallel::makePSOCKcluster(nslaves)
        doParallel::registerDoParallel(cl)
      }
    }
  }
  if("FP" %in% test){
    if(is.null(params$paramFP$basis)){ basis = "Fourier" }else{ basis = params$paramFP$basis }
    if(!(basis %in% c("Fourier", "b-spline", "own"))){ stop("argument params$paramFP$basis must be one of the following: 'Fourier', 'b-spline', 'own'") }
    if(is.null(params$paramFP$B.FP)){ nrFP = 1000 }else{ nrFP = params$paramFP$B.FP }
    if(nrFP < 1){ stop("invalid number of permutations (params$paramFP$B.FP)") }
    if(basis == "Fourier"){
      if(!("fda" %in% rownames(installed.packages()))){
        stop("Please install package 'fda'")
      }
      x = as.matrix(x); n = ncol(x); p = nrow(x)
      if(is.null(params$paramFP$minK)){ minK = 3 }else{ minK = params$paramFP$minK }
      if(minK < 1){ stop("invalid argument params$paramFP$minK") }
      if(is.null(params$paramFP$maxK)){  if(p %% 2 == 1){ maxK = p - 2 }else{ maxK = p - 1 }
      }else{ maxK = params$paramFP$maxK }
      if(maxK < 1){ stop("invalid argument params$paramFP$maxK") }
      if(p <= maxK){ if(p %% 2 == 1){ maxK = p - 2 }else{ maxK = p - 1 } }
      if((maxK %% 2) == 0){ stop("the maximum number of basis functions (params$paramFP$maxK) is even") }
      if(n != length(group.label)){
        stop("number of observations (number of columns in x) and number of elements
             in vector of group labels (group.label) must be the same")
      }
      if(is.null(params$paramFP$criterion)){ criterion = "BIC" }else{ criterion = params$paramFP$criterion }
      if(!(criterion %in% c("BIC", "eBIC", "AIC", "AICc", "NO"))){ stop("argument params$paramFP$criterion must be one of the following: 'BIC', 'eBIC', 'AIC', 'AICc', 'NO'") }
      if(criterion == "eBIC"){
        if(is.null(params$paramFP$gamma.eBIC)){ gamma.eBIC = 0.5 }else{ gamma.eBIC = params$paramFP$gamma.eBIC }
        if((gamma.eBIC < 0)|(gamma.eBIC > 1)){
          stop("argument params$paramFP$gamma.eBIC must belong to [0,1]")
        }
      }
      if(criterion != "NO"){
        if(is.null(params$paramFP$method)){ method = "mode" }else{ method = params$paramFP$method }
        if(!(method %in% c("mode", "min", "max", "mean"))){ stop("argument params$paramFP$method must be one of the following: 'mode', 'min', 'max', 'mean'") }
        if((minK %% 2) == 0){ stop("the minimum number of basis functions (params$paramFP$minK) is even") }
        if(is.null(params$paramFP$int)){
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(seq(minK, maxK, 2)))
            for(i in seq(minK, maxK, 2)){
              fbasis = fda::create.fourier.basis(c(0, p), i)
              fdata = fda::smooth.basisPar(1:p, x, fbasis)
              v[, (i+2-minK)/2] = switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*log(p),
                                         "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*(log(p) + 2*gamma.eBIC*log(maxK)),
                                         "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i,
                                         "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = seq(minK, maxK, 2), .combine = cbind) %dopar%
            {
              fbasis = fda::create.fourier.basis(c(0, p), i)
              fdata = fda::smooth.basisPar(1:p, x, fbasis)
              switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*log(p),
                     "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*(log(p) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i,
                     "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          KK = numeric(n)
          for(j in 1:n) KK[j] = 2*which(v[j,] == min(v[j,])) + minK - 2
          if(method == "mode"){
            temp = table(as.vector(KK))
            K = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(method == "min"){ K = min(KK) }
          if(method == "max"){ K = max(KK) }
          if(method == "mean"){ if(floor(mean(KK)) %% 2 == 1){ K = floor(mean(KK)) }else{ K = floor(mean(KK))-1 } }
        }else{
          if(length(params$paramFP$int) != 2){ stop("argument params$paramFP$int must be of length two") }
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(seq(minK, maxK, 2)))
            for(i in seq(minK, maxK, 2)){
              fbasis = fda::create.fourier.basis(rangeval = c(params$paramFP$int[1], params$paramFP$int[2]), nbasis = i)
              fdata = fda::smooth.basisPar(seq(params$paramFP$int[1], params$paramFP$int[2], length = p), x, fbasis)
              v[, (i+2-minK)/2] = switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*log(p),
                                         "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*(log(p) + 2*gamma.eBIC*log(maxK)),
                                         "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i,
                                         "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = seq(minK, maxK, 2), .combine = cbind) %dopar%
            {
              fbasis = fda::create.fourier.basis(rangeval = c(params$paramFP$int[1], params$paramFP$int[2]), nbasis = i)
              fdata = fda::smooth.basisPar(seq(params$paramFP$int[1], params$paramFP$int[2], length = p), x, fbasis)
              switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*log(p),
                     "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + i*(log(p) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i,
                     "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          KK = numeric(n)
          for(j in 1:n) KK[j] = 2*which(v[j,] == min(v[j,])) + minK - 2
          if(method == "mode"){
            temp = table(as.vector(KK))
            K = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(method == "min"){ K = min(KK) }
          if(method == "max"){ K = max(KK) }
          if(method == "mean"){ if(floor(mean(KK)) %% 2 == 1){ K = floor(mean(KK)) }else{ K = floor(mean(KK))-1 } }
        }
      }else{ K = maxK }

      if(is.null(params$paramFP$int)){
        fbasis = fda::create.fourier.basis(c(0, p), K)
        lfdata = vector("list", l)
        mfdata = matrix(0, nrow = K, ncol = n)
        for(i in 1:l){
          lfdata[[i]] = fda::Data2fd(1:p, x[, group.label == group.label0[i]], fbasis)$coefs
          mfdata[, group.label == group.label0[i]] = lfdata[[i]]
        }
      }else{
        if(length(params$paramFP$int) != 2){ stop("argument params$paramFP$int must be of length two") }
        fbasis = fda::create.fourier.basis(c(params$paramFP$int[1], params$paramFP$int[2]), K)
        lfdata = vector("list", l)
        mfdata = matrix(0, nrow = K, ncol = n)
        for(i in 1:l){
          lfdata[[i]] = fda::Data2fd(seq(params$paramFP$int[1], params$paramFP$int[2], length = p), x[, group.label == group.label0[i]], fbasis)$coefs
          mfdata[, group.label == group.label0[i]] = lfdata[[i]]
        }
      }

      a = 0; b = 0; c = 0
      for(i in 1:l){
        a = a + sum(t(lfdata[[i]]) %*% lfdata[[i]])/n.i[i]
        for(j in 1:l) b = b + sum(t(lfdata[[i]]) %*% (lfdata[[j]]))
        c = c + sum((lfdata[[i]])^2)
      }
      b = b/n; statFP0 = (a-b)/(c-a); statFP = statFP0*(n-l)/(l-1)
      if(parallel.method == "parallel.method0"){
        FP.p = numeric(nrFP)
        for(ii in 1:nrFP){
          a = numeric(l)
          Perm = sample(1:n)
          for(j in 1:l){
            if(j == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(j-1)])+1):sum(n.i[1:j])] }
            a[j] = sum(t(mfdata[, l1]) %*% mfdata[, l1])/n.i[j]
          }
          FP.p[ii] = (sum(a)-b)/(c-sum(a))
        }
      }
      if(parallel.method == "parallel.method1"){
        FP.p = foreach(ii = 1:nrFP, .combine = "c") %dopar%
        {
          a = numeric(l)
          Perm = sample(1:n)
          for(j in 1:l){
            if(j == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(j-1)])+1):sum(n.i[1:j])] }
            a[j] = sum(t(mfdata[, l1]) %*% mfdata[, l1])/n.i[j]
          }
          (sum(a)-b)/(c-sum(a))
        }
      }
      pvalueFP = mean(FP.p >= statFP0)
      }
    if(basis == "b-spline"){
      if(!("fda" %in% rownames(installed.packages()))){
        stop("Please install package 'fda'")
      }
      x = as.matrix(x); n = ncol(x); p = nrow(x)
      if(is.null(params$paramFP$norder)){ norder = 4 }else{ norder = params$paramFP$norder }
      if(norder < 1){ stop("invalid argument params$paramFP$norder") }
      if(is.null(params$paramFP$minK)){ minK = norder }else{ minK = params$paramFP$minK }
      if(minK < 1){ stop("invalid argument params$paramFP$minK") }
      if(minK < norder){ minK = norder }
      if(is.null(params$paramFP$maxK)){ maxK = p }else{ maxK = params$paramFP$maxK }
      if(maxK < 1){ stop("invalid argument params$paramFP$maxK") }
      if(p <= maxK){ maxK = p }
      if(n != length(group.label)){
        stop("number of observations (number of columns in x) and number of elements
             in vector of group labels (group.label) must be the same")
      }
      if(is.null(params$paramFP$criterion)){ criterion = "BIC" }else{ criterion = params$paramFP$criterion }
      if(!(criterion %in% c("BIC", "eBIC", "AIC", "AICc", "NO"))){ stop("argument params$paramFP$criterion must be one of the following: 'BIC', 'eBIC', 'AIC', 'AICc', 'NO'") }
      if(criterion == "eBIC"){
        if(is.null(params$paramFP$gamma.eBIC)){ gamma.eBIC = 0.5 }else{ gamma.eBIC = params$paramFP$gamma.eBIC }
        if((gamma.eBIC < 0)|(gamma.eBIC > 1)){
          stop("argument params$paramFP$gamma.eBIC must belong to [0,1]")
        }
      }
      if(criterion != "NO"){
        if(is.null(params$paramFP$method)){ method = "mode" }else{ method = params$paramFP$method }
        if(!(method %in% c("mode", "min", "max", "mean"))){ stop("argument params$paramFP$method must be one of the following: 'mode', 'min', 'max', 'mean'") }
        mmm = minK:maxK
        if(is.null(params$paramFP$int)){
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(mmm))
            for(i in 1:length(mmm)){
              fbasis = fda::create.bspline.basis(c(0, p), mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(1:p, x, fbasis)
              v[, i] = switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*log(p),
                              "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*(log(p) + 2*gamma.eBIC*log(maxK)),
                              "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i],
                              "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = 1:length(mmm), .combine = cbind) %dopar%
            {
              fbasis = fda::create.bspline.basis(c(0, p), mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(1:p, x, fbasis)
              switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*log(p),
                     "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*(log(p) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i],
                     "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          KK = numeric(n)
          for(j in 1:n) KK[j] = mmm[which(v[j,] == min(v[j,]))]
          if(method == "mode"){
            temp = table(as.vector(KK))
            K = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(method == "min"){ K = min(KK) }
          if(method == "max"){ K = max(KK) }
          if(method == "mean"){ K = floor(mean(KK)) }
        }else{
          if(length(params$paramFP$int) != 2){ stop("argument params$paramFP$int must be of length two") }
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(mmm))
            for(i in 1:length(mmm)){
              fbasis = fda::create.bspline.basis(rangeval = c(params$paramFP$int[1], params$paramFP$int[2]), nbasis = mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(seq(params$paramFP$int[1], params$paramFP$int[2], length = p), x, fbasis)
              v[, i] = switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*log(p),
                              "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*(log(p) + 2*gamma.eBIC*log(maxK)),
                              "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i],
                              "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = 1:length(mmm), .combine = cbind) %dopar%
            {
              fbasis = fda::create.bspline.basis(rangeval = c(params$paramFP$int[1], params$paramFP$int[2]), nbasis = mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(seq(params$paramFP$int[1], params$paramFP$int[2], length = p), x, fbasis)
              switch(criterion, "BIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*log(p),
                     "eBIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + mmm[i]*(log(p) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i],
                     "AICc" = p*log(p*(1-fdata$df/p)^2*fdata$gcv/p) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          KK = numeric(n)
          for(j in 1:n) KK[j] = mmm[which(v[j,] == min(v[j,]))]
          if(method == "mode"){
            temp = table(as.vector(KK))
            K = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(method == "min"){ K = min(KK) }
          if(method == "max"){ K = max(KK) }
          if(method == "mean"){ K = floor(mean(KK)) }
        }
      }else{ K = maxK }

      if(is.null(params$paramFP$int)){
        fbasis = fda::create.bspline.basis(c(0, p), K, norder = norder)
        lfdata = vector("list", l)
        mfdata = matrix(0, nrow = K, ncol = n)
        for(i in 1:l){
          lfdata[[i]] = fda::Data2fd(1:p, x[, group.label == group.label0[i]], fbasis)$coefs
          mfdata[, group.label == group.label0[i]] = lfdata[[i]]
        }
      }else{
        if(length(params$paramFP$int) != 2){ stop("argument params$paramFP$int must be of length two") }
        fbasis = fda::create.bspline.basis(c(params$paramFP$int[1], params$paramFP$int[2]), K, norder = norder)
        lfdata = vector("list", l)
        mfdata = matrix(0, nrow = K, ncol = n)
        for(i in 1:l){
          lfdata[[i]] = fda::Data2fd(seq(params$paramFP$int[1], params$paramFP$int[2], length = p), x[, group.label == group.label0[i]], fbasis)$coefs
          mfdata[, group.label == group.label0[i]] = lfdata[[i]]
        }
      }

      bspline.cross.prod.matrix = fda::inprod(fbasis, fbasis)
      a = 0; b = 0; c = 0
      for(i in 1:l){
        a = a + sum(t(lfdata[[i]]) %*% bspline.cross.prod.matrix %*% lfdata[[i]])/n.i[i]
        for(j in 1:l) b = b + sum(t(lfdata[[i]]) %*% bspline.cross.prod.matrix %*% (lfdata[[j]]))
        c = c + sum(diag(t(lfdata[[i]]) %*% bspline.cross.prod.matrix %*% lfdata[[i]]))
      }
      b = b/n; statFP0 = (a-b)/(c-a); statFP = statFP0*(n-l)/(l-1)
      if(parallel.method == "parallel.method0"){
        FP.p = numeric(nrFP)
        for(ii in 1:nrFP){
          a = numeric(l)
          Perm = sample(1:n)
          for(j in 1:l){
            if(j == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(j-1)])+1):sum(n.i[1:j])] }
            a[j] = sum(t(mfdata[, l1]) %*% bspline.cross.prod.matrix %*% mfdata[, l1])/n.i[j]
          }
          FP.p[ii] = (sum(a)-b)/(c-sum(a))
        }
      }
      if(parallel.method == "parallel.method1"){
        FP.p = foreach(ii = 1:nrFP, .combine = "c") %dopar%
        {
          a = numeric(l)
          Perm = sample(1:n)
          for(j in 1:l){
            if(j == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(j-1)])+1):sum(n.i[1:j])] }
            a[j] = sum(t(mfdata[, l1]) %*% bspline.cross.prod.matrix %*% mfdata[, l1])/n.i[j]
          }
          (sum(a)-b)/(c-sum(a))
        }
      }
      pvalueFP = mean(FP.p >= statFP0)
      }
    if(basis == "own"){
      own.basis = params$paramFP$own.basis
      mfdata = own.basis
      n = ncol(mfdata); K = nrow(mfdata)
      if(n != length(group.label)){ stop("the number of observations is not equal to the number of elements in vector of labels") }
      lfdata = vector("list", l)
      for(i in 1:l){ lfdata[[i]] = mfdata[, group.label == group.label0[i]] }

      own.cross.prod.mat = params$paramFP$own.cross.prod.mat
      if(K != nrow(own.cross.prod.mat)){ stop("invalid argument params$paramFP$own.cross.prod.mat") }
      if(K != ncol(own.cross.prod.mat)){ stop("invalid argument params$paramFP$own.cross.prod.mat") }

      a = 0; b = 0; c = 0
      for(i in 1:l){
        a = a + sum(t(lfdata[[i]]) %*% own.cross.prod.mat %*% lfdata[[i]])/n.i[i]
        for(j in 1:l) b = b + sum(t(lfdata[[i]]) %*% own.cross.prod.mat %*% (lfdata[[j]]))
        c = c + sum(diag(t(lfdata[[i]]) %*% own.cross.prod.mat %*% lfdata[[i]]))
      }
      b = b/n; statFP0 = (a-b)/(c-a); statFP = statFP0*(n-l)/(l-1)
      if(parallel.method == "parallel.method0"){
        FP.p = numeric(nrFP)
        for(ii in 1:nrFP){
          a = numeric(l)
          Perm = sample(1:n)
          for(j in 1:l){
            if(j == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(j-1)])+1):sum(n.i[1:j])] }
            a[j] = sum(t(mfdata[, l1]) %*% own.cross.prod.mat %*% mfdata[, l1])/n.i[j]
          }
          FP.p[ii] = (sum(a)-b)/(c-sum(a))
        }
      }
      if(parallel.method == "parallel.method1"){
        FP.p = foreach(ii = 1:nrFP, .combine = "c") %dopar%
        {
          a = numeric(l)
          Perm = sample(1:n)
          for(j in 1:l){
            if(j == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(j-1)])+1):sum(n.i[1:j])] }
            a[j] = sum(t(mfdata[, l1]) %*% own.cross.prod.mat %*% mfdata[, l1])/n.i[j]
          }
          (sum(a)-b)/(c-sum(a))
        }
      }
      pvalueFP = mean(FP.p >= statFP0)
    }
    if(basis == "own"){
      resultFP = list(statFP = statFP, pvalueFP = pvalueFP, B.FP = nrFP, basis = basis, own.basis = own.basis, own.cross.prod.mat = own.cross.prod.mat)
    }else{
      if(criterion == "NO"){
        if(is.null(params$paramFP$int)){
          if(basis == "Fourier"){
            resultFP = list(statFP = statFP, pvalueFP = pvalueFP, B.FP = nrFP, basis = basis, criterion = criterion, K = K, maxK = maxK)
          }
          if(basis == "b-spline"){
            resultFP = list(statFP = statFP, pvalueFP = pvalueFP, B.FP = nrFP, basis = basis, criterion = criterion, K = K, maxK = maxK, norder = norder)
          }
        }else{
          if(basis == "Fourier"){
            resultFP = list(statFP = statFP, pvalueFP = pvalueFP, int = params$paramFP$int, B.FP = nrFP, basis = basis, criterion = criterion, K = K, maxK = maxK)
          }
          if(basis == "b-spline"){
            resultFP = list(statFP = statFP, pvalueFP = pvalueFP, int = params$paramFP$int, B.FP = nrFP, basis = basis, criterion = criterion, K = K, maxK = maxK, norder = norder)
          }
        }
      }else{
        if(is.null(params$paramFP$int)){
          if(criterion == "eBIC"){
            if(basis == "Fourier"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK, gamma.eBIC = gamma.eBIC)
            }
            if(basis == "b-spline"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK, norder = norder, gamma.eBIC = gamma.eBIC)
            }
          }else{
            if(basis == "Fourier"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK)
            }
            if(basis == "b-spline"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK, norder = norder)
            }
          }
        }else{
          if(criterion == "eBIC"){
            if(basis == "Fourier"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, int = params$paramFP$int, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK, gamma.eBIC = gamma.eBIC)
            }
            if(basis == "b-spline"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, int = params$paramFP$int, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK, norder = norder, gamma.eBIC = gamma.eBIC)
            }
          }else{
            if(basis == "Fourier"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, int = params$paramFP$int, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK)
            }
            if(basis == "b-spline"){
              resultFP = list(statFP = statFP, pvalueFP = pvalueFP, int = params$paramFP$int, B.FP = nrFP, basis = basis, criterion = criterion, method = method, K = K, minK = minK, maxK = maxK, norder = norder)
            }
          }
        }
      }
    }
  }
  if(any(c("CH", "CS") %in% test)){
    if(!("MASS" %in% rownames(installed.packages()))){
      stop("Please install package 'MASS'")
    }
    x = as.matrix(x); n = ncol(x)
    if(n != length(group.label)){
      stop("number of observations (number of columns in x) and number of elements
           in vector of group labels (group.label) must be the same")
    }
    statCHCS = 0
    for(i in 1:(l-1)){
      for(j in (i+1):l){
        statCHCS = statCHCS + n.i[i]*sum((rowMeans(x[, group.label == group.label0[i]]) -
                                            rowMeans(x[, group.label == group.label0[j]]))^2)
      }
    }
    if("CH" %in% test){
      if(is.null(params$paramCH)){ nrCH = 10000 }else{ nrCH = params$paramCH }
      if(nrCH < 1){ stop("invalid number of discretized artificial trajectories for each process Z_i(t) (params$paramCH)") }
      cov.fun = ((n-1)/(n-l))*var(t(x)); Z = vector("list", l)
      for(i in 1:l) Z[[i]] = MASS::mvrnorm(nrCH, numeric(nrow(cov.fun)), cov.fun)
      if(parallel.method == "parallel.method0"){
        V = numeric(nrCH)
        for(m in 1:nrCH){
          for(i in 1:(l-1)){
            for(j in (i+1):l){
              V[m] = V[m] + sum((Z[[i]][m,]-sqrt(n.i[i]/n.i[j])*Z[[j]][m,])^2)
            }
          }
        }
      }
      if(parallel.method == "parallel.method1"){
        V = foreach(m = 1:nrCH, .combine = "c") %dopar%
        {
          V = 0
          for(i in 1:(l-1)){
            for(j in (i+1):l){
              V = V + sum((Z[[i]][m,]-sqrt(n.i[i]/n.i[j])*Z[[j]][m,])^2)
            }
          }
          V
        }
      }
      resultCH = list(statCH = statCHCS, pvalueCH = mean(V > statCHCS), paramCH = nrCH)
    }
    if("CS" %in% test){
      if(is.null(params$paramCS)){ nrCS = 10000 }else{ nrCS = params$paramCS }
      if(nrCS < 1){ stop("invalid number of discretized artificial trajectories for each process Z_i(t) (params$paramCS)") }
      Z = vector("list", l)
      for(i in 1:l){
        cov.fun = var(t(x[, group.label == group.label0[i]]))
        Z[[i]] = MASS::mvrnorm(nrCS, numeric(nrow(cov.fun)), cov.fun)
      }
      if(parallel.method == "parallel.method0"){
        V = numeric(nrCS)
        for(m in 1:nrCS){
          for(i in 1:(l-1)){
            for(j in (i+1):l){
              V[m] = V[m] + sum((Z[[i]][m,]-sqrt(n.i[i]/n.i[j])*Z[[j]][m,])^2)
            }
          }
        }
      }
      if(parallel.method == "parallel.method1"){
        V = foreach(m = 1:nrCS, .combine = "c") %dopar%
        {
          V = 0
          for(i in 1:(l-1)){
            for(j in (i+1):l){
              V = V + sum((Z[[i]][m,]-sqrt(n.i[i]/n.i[j])*Z[[j]][m,])^2)
            }
          }
          V
        }
      }
      resultCS = list(statCS = statCHCS, pvalueCS = mean(V > statCHCS), paramCS = nrCS)
    }
    }
  if(any(c("L2N", "L2B", "L2b", "FN", "FB", "Fb", "GPF", "Fmaxb") %in% test)){
    x = as.matrix(x); n = ncol(x); p = nrow(x)
    if(n != length(group.label)){
      stop("number of observations (number of columns in x) and number of elements
           in vector of group labels (group.label) must be the same")
    }
    mu0 = rowMeans(x)
    vmu = matrix(0, nrow = l, ncol = p)
    z = matrix(0, nrow = n, ncol = p)
    SSR = 0; SSE = 0
    for(i in 1:l){
      xi = x[, group.label == group.label0[i]]
      mui = rowMeans(xi); vmu[i,] = mui
      zi = t(xi) - as.matrix(rep(1, n.i[i])) %*% mui
      if(i==1){ z[1:n.i[i],] = zi }else{ z[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = zi }
      SSR = SSR + n.i[i]*(mui - mu0)^2
      SSE = SSE + colSums(zi^2)
    }
    if(n > p){ gamma.z = t(z) %*% z/(n-l) }else{ gamma.z = z %*% t(z)/(n-l) }
    A = sum(diag(gamma.z)); B = sum(diag(gamma.z %*% gamma.z))
    A2N = A^2; B2N = B
    A2B = (n-l)*(n-l+1)/(n-l-1)/(n-l+2)*(A^2-2*B/(n-l+1))
    B2B = (n-l)^2/(n-l-1)/(n-l+2)*(B-A^2/(n-l))

    if(any(c("L2N", "L2B", "L2b") %in% test)){
      statL2 = sum(SSR)
      if("L2N" %in% test){
        betaL2N = B2N/A; kappaL2N = A2N/B2N
        pvalueL2N = 1-pchisq(statL2/betaL2N, (l-1)*kappaL2N)
        resultL2N = list(statL2 = statL2, pvalueL2N = pvalueL2N, betaL2N = betaL2N, dL2N = (l-1)*kappaL2N)
      }
      if("L2B" %in% test){
        betaL2B = B2B/A; kappaL2B = A2B/B2B
        pvalueL2B = 1-pchisq(statL2/betaL2B, (l-1)*kappaL2B)
        resultL2B = list(statL2 = statL2, pvalueL2B = pvalueL2B, betaL2B = betaL2B, dL2B = (l-1)*kappaL2B)
      }
      if("L2b" %in% test){
        if(is.null(params$paramL2b)){ nrL2b = 10000 }else{ nrL2b = params$paramL2b }
        if(nrL2b < 1){ stop("invalid number of bootstrap replicates (params$paramL2b)") }
        if(parallel.method == "parallel.method0"){
          statL2boot = numeric(nrL2b)
          for(ii in 1:nrL2b){
            vmuboot = matrix(0, nrow = l, ncol = p)
            for(i in 1:l){
              xi = x[, group.label == group.label0[i]]
              xiboot = xi[, floor(runif(n.i[i])*(n.i[i]-1)) + 1]
              vmuboot[i,] = rowMeans(xiboot)-vmu[i,]
            }
            mu0boot = n.i %*% vmuboot/n
            SSRboot = 0
            for(i in 1:l) SSRboot = SSRboot + n.i[i]*(vmuboot[i,]-mu0boot)^2
            statL2boot[ii] = sum(SSRboot)
          }
        }
        if(parallel.method == "parallel.method1"){
          statL2boot = foreach(ii = 1:nrL2b, .combine = "c") %dopar%
          {
            vmuboot = matrix(0, nrow = l, ncol = p)
            for(i in 1:l){
              xi = x[, group.label == group.label0[i]]
              xiboot = xi[, floor(runif(n.i[i])*(n.i[i]-1)) + 1]
              vmuboot[i,] = rowMeans(xiboot)-vmu[i,]
            }
            mu0boot = n.i %*% vmuboot/n
            SSRboot = 0
            for(i in 1:l) SSRboot = SSRboot + n.i[i]*(vmuboot[i,]-mu0boot)^2
            sum(SSRboot)
          }
        }
        pvalueL2b = mean(statL2boot >= statL2)
        resultL2b = list(statL2 = statL2, pvalueL2b = pvalueL2b, paramL2b = nrL2b)
      }
    }
    if(any(c("FN", "FB", "Fb") %in% test)){
      statF = sum(SSR)/sum(SSE)*(n-l)/(l-1)
      if("FN" %in% test){
        kappaFN = A2N/B2N
        pvalueFN = 1-pf(statF, (l-1)*kappaFN, (n-l)*kappaFN)
        resultFN = list(statF = statF, pvalueFN = pvalueFN, d1FN = (l-1)*kappaFN, d2FN = (n-l)*kappaFN)
      }
      if("FB" %in% test){
        kappaFB = A2B/B2B
        pvalueFB = 1-pf(statF, (l-1)*kappaFB, (n-l)*kappaFB)
        resultFB = list(statF = statF, pvalueFB = pvalueFB, d1FB = (l-1)*kappaFB, d2FB = (n-l)*kappaFB)
      }
      if("Fb" %in% test){
        if(is.null(params$paramFb)){ nrFb = 10000 }else{ nrFb = params$paramFb }
        if(nrFb < 1){ stop("invalid number of bootstrap replicates (params$paramFb)") }
        if(parallel.method == "parallel.method0"){
          statFboot = numeric(nrFb)
          for(ii in 1:nrFb){
            vmuboot = matrix(0, nrow = l, ncol = p)
            zboot = matrix(0, nrow = n, ncol = p)
            for(i in 1:l){
              xi = x[, group.label == group.label0[i]]
              xiboot = xi[, floor(runif(n.i[i])*(n.i[i]-1)) + 1]
              muiboot = rowMeans(xiboot)
              ziboot = t(xiboot) - as.matrix(rep(1, n.i[i])) %*% muiboot
              if(i==1){ zboot[1:n.i[i],] = ziboot }else{ zboot[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = ziboot }
              vmuboot[i,] = rowMeans(xiboot)-vmu[i,]
            }
            mu0boot = n.i %*% vmuboot/n
            SSRboot = 0
            for(i in 1:l) SSRboot = SSRboot + n.i[i]*(vmuboot[i,]-mu0boot)^2
            if(n > p){ gammaboot = t(zboot) %*% zboot/(n-l) }else{ gammaboot = zboot %*% t(zboot)/(n-l) }
            statFboot[ii] = sum(SSRboot)/sum(diag(gammaboot))/(l-1)
          }
        }
        if(parallel.method == "parallel.method1"){
          statFboot = foreach(ii = 1:nrFb, .combine = "c") %dopar%
          {
            vmuboot = matrix(0, nrow = l, ncol = p)
            zboot = matrix(0, nrow = n, ncol = p)
            for(i in 1:l){
              xi = x[, group.label == group.label0[i]]
              xiboot = xi[, floor(runif(n.i[i])*(n.i[i]-1)) + 1]
              muiboot = rowMeans(xiboot)
              ziboot = t(xiboot) - as.matrix(rep(1, n.i[i])) %*% muiboot
              if(i==1){ zboot[1:n.i[i],] = ziboot }else{ zboot[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = ziboot }
              vmuboot[i,] = rowMeans(xiboot)-vmu[i,]
            }
            mu0boot = n.i %*% vmuboot/n
            SSRboot = 0
            for(i in 1:l) SSRboot = SSRboot + n.i[i]*(vmuboot[i,]-mu0boot)^2
            if(n > p){ gammaboot = t(zboot) %*% zboot/(n-l) }else{ gammaboot = zboot %*% t(zboot)/(n-l) }
            sum(SSRboot)/sum(diag(gammaboot))/(l-1)
          }
        }
        pvalueFb = mean(statFboot >= statF)
        resultFb = list(statF = statF, pvalueFb = pvalueFb, paramFb = nrFb)
      }
    }
    if("GPF" %in% test){
      statGPF = mean(SSR/SSE*(n-l)/(l-1))
      z = z/(as.matrix(rep(1, n)) %*% sqrt(colSums(z^2)))
      if(n >= p){ z = t(z) %*% z }else{ z = z %*% t(z) }
      betaGPF = (sum(diag(z %*% z))/p^2/(l-1))/((n-l)/(n-l-2))
      dGPF = ((n-l)/(n-l-2))^2/(sum(diag(z %*% z))/p^2/(l-1))
      pvalueGPF = 1-pchisq(statGPF/betaGPF,dGPF)
      resultGPF = list(statGPF = statGPF, pvalueGPF = pvalueGPF, betaGPF = betaGPF, dGPF = dGPF)
    }
    if("Fmaxb" %in% test){
      if(is.null(params$paramFmaxb)){ nrFmaxb = 10000 }else{ nrFmaxb = params$paramFmaxb }
      if(nrFmaxb < 1){ stop("invalid number of bootstrap replicates (params$paramFmaxb)") }
      statFmax = max(SSR/SSE*(n-l)/(l-1))
      xs = matrix(0, nrow = n, ncol = p)
      for(i in 1:l){ xs[group.label == group.label0[i],] = t(x[, group.label == group.label0[i]]) - as.matrix(rep(1, n.i[i])) %*% vmu[i,] }
      if(parallel.method == "parallel.method0"){
        statFmaxboot = numeric(nrFmaxb)
        for(ii in 1:nrFmaxb){
          vmuboot = matrix(0, nrow = l, ncol = p)
          zboot = matrix(0, nrow = n, ncol = p)
          xboot = xs[sample(1:n, replace = TRUE),]
          mu0boot = colMeans(xboot)
          SSRboot = 0
          for(i in 1:l){
            xiboot = xboot[group.label == group.label0[i],]
            muiboot = colMeans(xiboot)
            vmuboot[i,] = muiboot
            ziboot = xiboot - as.matrix(rep(1, n.i[i])) %*% muiboot
            if(i==1){ zboot[1:n.i[i],] = ziboot }else{ zboot[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = ziboot }
            SSRboot = SSRboot + n.i[i]*(muiboot - mu0boot)^2
          }
          SSEboot = diag(t(zboot) %*% zboot)
          statFmaxboot[ii] = max(SSRboot/SSEboot*(n-l)/(l-1))
        }
      }
      if(parallel.method == "parallel.method1"){
        statFmaxboot = foreach(ii = 1:nrFmaxb, .combine = "c") %dopar%
        {
          vmuboot = matrix(0, nrow = l, ncol = p)
          zboot = matrix(0, nrow = n, ncol = p)
          xboot = xs[sample(1:n, replace = TRUE),]
          mu0boot = colMeans(xboot)
          SSRboot = 0
          for(i in 1:l){
            xiboot = xboot[group.label == group.label0[i],]
            muiboot = colMeans(xiboot)
            vmuboot[i,] = muiboot
            ziboot = xiboot - as.matrix(rep(1, n.i[i])) %*% muiboot
            if(i==1){ zboot[1:n.i[i],] = ziboot }else{ zboot[(cumsum(n.i)[i-1]+1):cumsum(n.i)[i],] = ziboot }
            SSRboot = SSRboot + n.i[i]*(muiboot - mu0boot)^2
          }
          SSEboot = diag(t(zboot) %*% zboot)
          max(SSRboot/SSEboot*(n-l)/(l-1))
        }
      }
      pvalueFmaxb = mean(statFmaxboot >= statFmax)
      resultFmaxb = list(statFmax = statFmax, pvalueFmaxb = pvalueFmaxb, paramFmaxb = nrFmaxb)
    }
    }
  if("TRP" %in% test){
    if(!("doBy" %in% rownames(installed.packages()))){
      stop("Please install package 'doBY'")
    }
    if(!("MASS" %in% rownames(installed.packages()))){
      stop("Please install package 'MASS'")
    }
    if(!("magic" %in% rownames(installed.packages()))){
      stop("Please install package 'magic'")
    }
    x = as.matrix(x); n = ncol(x); p = nrow(x)
    if(n != length(group.label)){
      stop("number of observations (number of columns in x) and number of elements
           in vector of group labels (group.label) must be the same")
    }
    if(is.null(params$paramTRP$projection)){ projection = "GAUSS" }else{ projection = params$paramTRP$projection }
    if(!(projection %in% c("GAUSS", "BM"))){ stop("argument params$paramTRP$projection must be one of the following: 'GAUSS', 'BM'") }
    if(is.null(params$paramTRP$permutation)){ permutationTRP = FALSE }else{ permutationTRP = params$paramTRP$permutation }
    if(! is.logical(permutationTRP)){ stop("argument permutation is not logical (params$paramTRP$permutation)") }
    if(is.null(params$paramTRP$B.TRP)){ nrTRP = 10000 }else{ nrTRP = params$paramTRP$B.TRP }
    if(nrTRP < 1){ stop("invalid number of permutations (params$paramTRP$B.TRP)") }
    if(is.null(params$paramTRP$k)){ k = 30 }else{ k = params$paramTRP$k }
    if(any(k < 1)){ stop("invalid number of projections (k)") }
    pvalues.anova = numeric(length(k))
    pvalues.ATS.simple = numeric(length(k))
    pvalues.WTPS = numeric(length(k))
    all.data.proj = list()
    iik = 0
    modulo = function(z){ sqrt(sum(z^2)) }
    for(ik in k){
      ik.data.proj = matrix(0, nrow = n, ncol = ik)
      if(projection == "GAUSS"){
        z = matrix(rnorm(p * ik), nrow = ik, ncol = p)
        modu = apply(z, 1, modulo)
        z = z/modu
        if(parallel.method == "parallel.method0"){
          anova.p = numeric(ik)
          ATS.simple.p = numeric(ik)
          WTPS.p = numeric(ik)
          for(j in 1:ik){
            ik.data.proj[, j] = t(x) %*% z[j,]
            data.proj = ik.data.proj[, j]
            perm.WTPS = matrix(0, nrow = n, ncol = nrTRP)
            for(i.perm in 1:nrTRP){
              perm.WTPS[, i.perm] = sample(1:n)
            }
            WTPS.p[j] = WTPSp(data.proj, group.label, group.label0, n, n.i, l, perm.WTPS = perm.WTPS, nrTRP)
            if(permutationTRP == FALSE){
              anova.p[j] = 1-pf(anova.statistic.quick(data.proj, group.label), l-1, n-l)
              ATS.simple.p[j] = ATS.simple(data.proj, group.label)[2]
            }else{
              anova.s = anova.statistic.quick(data.proj, group.label)
              ATS.simple.s = ATS.simple(data.proj, group.label)[1]
              anova.perm = numeric(nrTRP)
              ATS.simple.perm = numeric(nrTRP)
              for(i.perm in 1:nrTRP){
                anova.perm[i.perm] = anova.statistic.quick(data.proj, sample(group.label))
                ATS.simple.perm[i.perm] = ATS.simple(data.proj, sample(group.label))[1]
              }
              anova.p[j] = mean(anova.perm >= anova.s)
              ATS.simple.p[j] = mean(ATS.simple.perm >= ATS.simple.s)
            }
          }
          iik = iik + 1
          all.data.proj[[iik]] = ik.data.proj
          pvalues.anova[iik] = min(ik*anova.p[order(anova.p)]/1:ik)
          pvalues.ATS.simple[iik] = min(ik*ATS.simple.p[order(ATS.simple.p)]/1:ik)
          pvalues.WTPS[iik] = min(ik*WTPS.p[order(WTPS.p)]/1:ik)
        }
        if(parallel.method == "parallel.method1"){
          rs = foreach(j = 1:ik, .combine = rbind) %dopar%
          {
            ik.data.proj[, j] = t(x) %*% z[j,]
            data.proj = ik.data.proj[, j]
            perm.WTPS = matrix(0, nrow = n, ncol = nrTRP)
            for(i.perm in 1:nrTRP){
              perm.WTPS[, i.perm] = sample(1:n)
            }
            WTPS.p = WTPSp(data.proj, group.label, group.label0, n, n.i, l, perm.WTPS = perm.WTPS, nrTRP)
            if(permutationTRP == FALSE){
              c(1-pf(anova.statistic.quick(data.proj, group.label), l-1, n-l),
                ATS.simple(data.proj, group.label)[2], WTPS.p, data.proj)
            }else{
              anova.s = anova.statistic.quick(data.proj, group.label)
              ATS.simple.s = ATS.simple(data.proj, group.label)[1]
              anova.perm = numeric(nrTRP)
              ATS.simple.perm = numeric(nrTRP)
              for(i.perm in 1:nrTRP){
                anova.perm[i.perm] = anova.statistic.quick(data.proj, sample(group.label))
                ATS.simple.perm[i.perm] = ATS.simple(data.proj, sample(group.label))[1]
              }
              c(mean(anova.perm >= anova.s),
                mean(ATS.simple.perm >= ATS.simple.s), WTPS.p, data.proj)
            }
          }
          if(ik == 1) rs = matrix(rs, nrow = 1, ncol = n+3)
          iik = iik + 1
          if(ik == 1){
            all.data.proj[[iik]] = matrix(rs[, 4:(n+3)], ncol = 1)
          }else{
            all.data.proj[[iik]] = t(rs[, 4:(n+3)])
          }
          pvalues.anova[iik] = min(ik*(rs[, 1])[order(rs[, 1])]/1:ik)
          pvalues.ATS.simple[iik] = min(ik*(rs[, 2])[order(rs[, 2])]/1:ik)
          pvalues.WTPS[iik] = min(ik*(rs[, 3])[order(rs[, 3])]/1:ik)
        }
      }else{
        if(parallel.method == "parallel.method0"){
          anova.p = numeric(ik)
          ATS.simple.p = numeric(ik)
          WTPS.p = numeric(ik)
          for(j in 1:ik){
            bm.p = cumsum(c(mean(x[1,]), rnorm(p-1, mean = 0, sd = 1)))/sqrt(p)
            bm.p = bm.p/modulo(bm.p)
            ik.data.proj[, j] = t(x) %*% as.matrix(bm.p)
            data.proj = ik.data.proj[, j]
            perm.WTPS = matrix(0, nrow = n, ncol = nrTRP)
            for(i.perm in 1:nrTRP){
              perm.WTPS[, i.perm] = sample(1:n)
            }
            WTPS.p[j] = WTPSp(data.proj, group.label, group.label0, n, n.i, l, perm.WTPS = perm.WTPS, nrTRP)
            if(permutationTRP == FALSE){
              anova.p[j] = 1-pf(anova.statistic.quick(data.proj, group.label), l-1, n-l)
              ATS.simple.p[j] = ATS.simple(data.proj, group.label)[2]
            }else{
              anova.s = anova.statistic.quick(data.proj, group.label)
              ATS.simple.s = ATS.simple(data.proj, group.label)[1]
              anova.perm = numeric(nrTRP)
              ATS.simple.perm = numeric(nrTRP)
              for(i.perm in 1:nrTRP){
                anova.perm[i.perm] = anova.statistic.quick(data.proj, sample(group.label))
                ATS.simple.perm[i.perm] = ATS.simple(data.proj, sample(group.label))[1]
              }
              anova.p[j] = mean(anova.perm >= anova.s)
              ATS.simple.p[j] = mean(ATS.simple.perm >= ATS.simple.s)
            }
          }
          iik = iik + 1
          all.data.proj[[iik]] = ik.data.proj
          pvalues.anova[iik] = min(ik*anova.p[order(anova.p)]/1:ik)
          pvalues.ATS.simple[iik] = min(ik*ATS.simple.p[order(ATS.simple.p)]/1:ik)
          pvalues.WTPS[iik] = min(ik*WTPS.p[order(WTPS.p)]/1:ik)
        }
        if(parallel.method == "parallel.method1"){
          rs = foreach(j = 1:ik, .combine = rbind) %dopar%
          {
            bm.p = cumsum(c(mean(x[1,]), rnorm(p-1, mean = 0, sd = 1)))/sqrt(p)
            bm.p = bm.p/modulo(bm.p)
            ik.data.proj[, j] = t(x) %*% as.matrix(bm.p)
            data.proj = ik.data.proj[, j]
            perm.WTPS = matrix(0, nrow = n, ncol = nrTRP)
            for(i.perm in 1:nrTRP){
              perm.WTPS[, i.perm] = sample(1:n)
            }
            WTPS.p = WTPSp(data.proj, group.label, group.label0, n, n.i, l, perm.WTPS = perm.WTPS, nrTRP)
            if(permutationTRP == FALSE){
              c(1-pf(anova.statistic.quick(data.proj, group.label), l-1, n-l),
                ATS.simple(data.proj, group.label)[2], WTPS.p, data.proj)
            }else{
              anova.s = anova.statistic.quick(data.proj, group.label)
              ATS.simple.s = ATS.simple(data.proj, group.label)[1]
              anova.perm = numeric(nrTRP)
              ATS.simple.perm = numeric(nrTRP)
              for(i.perm in 1:nrTRP){
                anova.perm[i.perm] = anova.statistic.quick(data.proj, sample(group.label))
                ATS.simple.perm[i.perm] = ATS.simple(data.proj, sample(group.label))[1]
              }
              c(mean(anova.perm >= anova.s),
                mean(ATS.simple.perm >= ATS.simple.s), WTPS.p, data.proj)
            }
          }
          if(ik == 1) rs = matrix(rs, nrow = 1, ncol = n+3)
          iik = iik + 1
          if(ik == 1){
            all.data.proj[[iik]] = matrix(rs[, 4:(n+3)], ncol = 1)
          }else{
            all.data.proj[[iik]] = t(rs[, 4:(n+3)])
          }
          pvalues.anova[iik] = min(ik*(rs[, 1])[order(rs[, 1])]/1:ik)
          pvalues.ATS.simple[iik] = min(ik*(rs[, 2])[order(rs[, 2])]/1:ik)
          pvalues.WTPS[iik] = min(ik*(rs[, 3])[order(rs[, 3])]/1:ik)
        }
      }
    }
    resultTRP = list(pvalues.anova = pvalues.anova, pvalues.ATS = pvalues.ATS.simple,
                     pvalues.WTPS = pvalues.WTPS, data.projections = all.data.proj,
                     k = k, projection = projection,
                     permutation = permutationTRP, B.TRP = nrTRP)
    }
  if(any(c("FP", "CH", "CS", "L2b", "Fb", "Fmaxb", "TRP") %in% test)){
    if(parallel.method == "parallel.method1"){
      parallel::stopCluster(cl)
    }
  }

  list.results = list()
  for(ii in test){
    if(ii == "FP"){
      list.results[[which(test == ii)]] = resultFP
      names(list.results)[which(test == ii)] = "FP"
    }
    if(ii == "CH"){
      list.results[[which(test == ii)]] = resultCH
      names(list.results)[which(test == ii)] = "CH"
    }
    if(ii == "CS"){
      list.results[[which(test == ii)]] = resultCS
      names(list.results)[which(test == ii)] = "CS"
    }
    if(ii == "L2N"){
      list.results[[which(test == ii)]] = resultL2N
      names(list.results)[which(test == ii)] = "L2N"
    }
    if(ii == "L2B"){
      list.results[[which(test == ii)]] = resultL2B
      names(list.results)[which(test == ii)] = "L2B"
    }
    if(ii == "L2b"){
      list.results[[which(test == ii)]] = resultL2b
      names(list.results)[which(test == ii)] = "L2b"
    }
    if(ii == "FN"){
      list.results[[which(test == ii)]] = resultFN
      names(list.results)[which(test == ii)] = "FN"
    }
    if(ii == "FB"){
      list.results[[which(test == ii)]] = resultFB
      names(list.results)[which(test == ii)] = "FB"
    }
    if(ii == "Fb"){
      list.results[[which(test == ii)]] = resultFb
      names(list.results)[which(test == ii)] = "Fb"
    }
    if(ii == "GPF"){
      list.results[[which(test == ii)]] = resultGPF
      names(list.results)[which(test == ii)] = "GPF"
    }
    if(ii == "Fmaxb"){
      list.results[[which(test == ii)]] = resultFmaxb
      names(list.results)[which(test == ii)] = "Fmaxb"
    }
    if(ii == "TRP"){
      list.results[[which(test == ii)]] = resultTRP
      names(list.results)[which(test == ii)] = "TRP"
    }
  }
  list.results$data = x
  list.results$group.label = group.label
  list.results$parallel = parallel
  list.results$nslaves = nslaves
  class(list.results) = "fanovatests"
  return(list.results)
}
