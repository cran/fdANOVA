fmanova.ptbfr = function(x = NULL, group.label, int, B = 1000,
                         parallel = FALSE, nslaves = NULL,
                         basis = c("Fourier", "b-spline", "own"),
                         own.basis, own.cross.prod.mat,
                         criterion = c("BIC", "eBIC", "AIC", "AICc", "NO"),
                         commonK = c("mode", "min", "max", "mean"),
                         minK = NULL, maxK = NULL, norder = 4, gamma.eBIC = 0.5){
  group.label0 = unique(group.label)
  l = length(group.label0)
  n.i = numeric(l)
  for(i in 1:l) n.i[i] = sum(group.label == group.label0[i])
  if(! is.logical(parallel)){ stop("argument parallel is not logical") }
  basis = match.arg(basis)
  if(B < 1){ stop("invalid number of permutations (B)") }
  criterion = match.arg(criterion)
  commonK = match.arg(commonK)
  if(any(is.na(group.label))){ stop("argument group.label can not contain NA values") }
  if((gamma.eBIC < 0)|(gamma.eBIC > 1)){ stop("argument gamma.eBIC must belong to [0,1]") }
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
  if(basis == "Fourier"){
    if(!("fda" %in% rownames(installed.packages()))){
      stop("Please install package 'fda'")
    }
    n = ncol(x[[1]])
    p = length(x)
    I = nrow(x[[1]])
    if(n != length(group.label)){ stop("the number of observations is not equal to the number of elements in vector of labels") }
    if(is.null(minK)){ minK = 3 }
    if(minK < 1){ stop("invalid argument minK") }
    if(is.null(maxK)){  if(I %% 2 == 1){ maxK = I - 2 }else{ maxK = I - 1 } }
    if(maxK < 1){ stop("invalid argument maxK") }
    if(I <= maxK){ if(I %% 2 == 1) maxK = I - 2 else maxK = I - 1 }
    if((maxK %% 2) == 0){ stop("the maximum number of basis functions (maxK) is even") }
    if(criterion != "NO"){
      if((minK %% 2) == 0){ stop("the minimum number of basis functions (minK) is even") }
      KK = numeric(p)
      if(missing(int)){
        for(j in 1:p){
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(seq(minK, maxK, 2)))
            for(i in seq(minK, maxK, 2)){
              fbasis = fda::create.fourier.basis(c(0, I), i)
              fdata = fda::smooth.basisPar(1:I, x[[j]], fbasis)
              v[, (i+2-minK)/2] = switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*log(I),
                                         "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*(log(I) + 2*gamma.eBIC*log(maxK)),
                                         "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i,
                                         "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = seq(minK, maxK, 2), .combine = cbind) %dopar%
            {
              fbasis = fda::create.fourier.basis(c(0, I), i)
              fdata = fda::smooth.basisPar(1:I, x[[j]], fbasis)
              switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*log(I),
                     "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*(log(I) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i,
                     "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          K = numeric(n)
          for(jj in 1:n) K[jj] = 2*which(v[jj,] == min(v[jj,])) + minK - 2
          if(commonK == "mode"){
            temp = table(as.vector(K))
            KK[j] = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(commonK == "min"){ KK[j] = min(K) }
          if(commonK == "max"){ KK[j] = max(K) }
          if(commonK == "mean"){ if(floor(mean(K)) %% 2 == 1){ KK[j] = floor(mean(K)) }else{ KK[j] = floor(mean(K))-1 } }
        }
      }else{
        if(length(int) != 2){ stop("argument int must be of length two") }
        for(j in 1:p){
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(seq(minK, maxK, 2)))
            for(i in seq(minK, maxK, 2)){
              fbasis = fda::create.fourier.basis(rangeval = c(int[1], int[2]), nbasis = i)
              fdata = fda::smooth.basisPar(seq(int[1], int[2], length = I), x[[j]], fbasis)
              v[, (i+2-minK)/2] = switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*log(I),
                                         "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*(log(I) + 2*gamma.eBIC*log(maxK)),
                                         "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i,
                                         "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = seq(minK, maxK, 2), .combine = cbind) %dopar%
            {
              fbasis = fda::create.fourier.basis(rangeval = c(int[1], int[2]), nbasis = i)
              fdata = fda::smooth.basisPar(seq(int[1], int[2], length = I), x[[j]], fbasis)
              switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*log(I),
                     "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + i*(log(I) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i,
                     "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*i + (2 * i * (i + 1))/(n - i - 1))
            }
          }
          K = numeric(n)
          for(jj in 1:n) K[jj] = 2*which(v[jj,] == min(v[jj,])) + minK - 2
          if(commonK == "mode"){
            temp = table(as.vector(K))
            KK[j] = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(commonK == "min"){ KK[j] = min(K) }
          if(commonK == "max"){ KK[j] = max(K) }
          if(commonK == "mean"){ if(floor(mean(K)) %% 2 == 1){ KK[j] = floor(mean(K)) }else{ KK[j] = floor(mean(K))-1 } }
        }
      }
    }else{ KK = rep(maxK, p) }
    KM = max(KK)

    data.fd = vector("list", p)
    if(missing(int)){
      for(i in 1:p){
        fbasis = fda::create.fourier.basis(c(0, I), KK[i])
        data.fd[[i]] = t(fda::Data2fd(1:I, x[[i]], fbasis)$coefs)
      }
    }else{
      if(length(int) != 2){ stop("argument int must be of length two") }
      for(i in 1:p){
        fbasis = fda::create.fourier.basis(c(int[1], int[2]), KK[i])
        data.fd[[i]] = t(fda::Data2fd(seq(int[1], int[2], length = I), x[[i]], fbasis)$coefs)
      }
    }

    data.ALFA = vector("list", n)
    for(i in 1:l){
      temp = vector("list", n.i[i])
      for(j in 1:n.i[i]){
        temp[[j]] = matrix(rep(0, p*KM), nrow = p)
        for(h in 1:p) temp[[j]][h, 1:KK[h]] =
            data.fd[[h]][group.label == group.label0[i],][j,]
      }
      data.ALFA[group.label == group.label0[i]] = temp
    }

    AA = 0; BB = 0; CC = 0
    for(i in 1:l){
      z1 = data.ALFA[group.label == group.label0[i]]
      z2 = lapply(z1, t)
      BB1 = 0
      for(j1 in 1:n.i[i]){
        AA = AA + z1[[j1]] %*% z2[[j1]]
        for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% z2[[j2]]
        for(m in 1:l){
          z3 = lapply(data.ALFA[group.label == group.label0[m]], t)
          for(j3 in 1:n.i[m]) CC = CC + z1[[j1]] %*% z3[[j3]]
        }
      }
      BB = BB + BB1/n.i[i]
    }
    CC = CC/n
    Y = solve(AA - CC)
    W.0 = det(AA - BB)
    LH.0 = sum(diag((BB - CC) %*% solve(AA - BB)))
    P.0 = sum(diag((BB - CC) %*% Y))
    R.0 = max(Re(eigen((BB - CC) %*% solve(AA - BB))$values))

    W.P = numeric(B); LH.P = numeric(B); P.P = numeric(B); R.P = numeric(B)
    if(parallel.method == "parallel.method0"){
      for(jP in 1:B){
        BB = 0
        Perm = sample(1:n)
        for(i in 1:l){
          if(i == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(i-1)])+1):sum(n.i[1:i])] }
          z1 = data.ALFA[l1]
          z2 = lapply(z1, t)
          BB1 = 0
          for(j1 in 1:n.i[i]) for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% z2[[j2]]
          BB = BB + BB1/n.i[i]
        }
        W.P[jP] = det(AA - BB)
        LH.P[jP] = sum(diag((BB - CC) %*% solve(AA - BB)))
        P.P[jP] = sum(diag((BB - CC) %*% Y))
        R.P[jP] = max(Re(eigen((BB - CC) %*% solve(AA - BB))$values))
      }
      W = W.0/det(AA - CC); pvalue.W = mean(W.P <= W.0); LH = LH.0; pvalue.LH = mean(LH.P >= LH.0)
      P = P.0; pvalue.P = mean(P.P >= P.0); R = R.0; pvalue.R = mean(R.P >= R.0)
    }
    if(parallel.method == "parallel.method1"){
      rs.perm = foreach(jP = 1:B, .combine = rbind) %dopar%
      {
        BB = 0
        Perm = sample(1:n)
        for(i in 1:l){
          if(i == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(i-1)])+1):sum(n.i[1:i])] }
          z1 = data.ALFA[l1]
          z2 = lapply(z1, t)
          BB1 = 0
          for(j1 in 1:n.i[i]) for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% z2[[j2]]
          BB = BB + BB1/n.i[i]
        }
        c(det(AA-BB), sum(diag((BB - CC) %*% solve(AA - BB))),
          sum(diag((BB - CC) %*% Y)), max(Re(eigen((BB - CC) %*% solve(AA - BB))$values)))
      }
      W = W.0/det(AA - CC); pvalue.W = mean(rs.perm[,1] <= W.0); LH = LH.0; pvalue.LH = mean(rs.perm[,2] >= LH.0)
      P = P.0; pvalue.P = mean(rs.perm[,3] >= P.0); R = R.0; pvalue.R = mean(rs.perm[,4] >= R.0)
      parallel::stopCluster(cl)
    }
  }
  if(basis == "b-spline"){
    if(!("fda" %in% rownames(installed.packages()))){
      stop("Please install package 'fda'")
    }
    n = ncol(x[[1]])
    p = length(x)
    I = nrow(x[[1]])
    if(n != length(group.label)){ stop("the number of observations is not equal to the number of elements in vector of labels") }
    if(norder < 1){ stop("invalid argument norder") }
    if(is.null(minK)){ minK = norder }
    if(minK < 1){ stop("invalid argument minK") }
    if(minK < norder){ minK = norder }
    if(is.null(maxK)){ maxK = I }
    if(maxK < 1){ stop("invalid argument maxK") }
    if(I <= maxK){ maxK = I }
    if(criterion != "NO"){
      KK = numeric(p)
      mmm = minK:maxK
      if(missing(int)){
        for(j in 1:p){
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(mmm))
            for(i in 1:length(mmm)){
              fbasis = fda::create.bspline.basis(c(0, I), mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(1:I, x[[j]], fbasis)
              v[, i] = switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*log(I),
                              "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*(log(I) + 2*gamma.eBIC*log(maxK)),
                              "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i],
                              "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = 1:length(mmm), .combine = cbind) %dopar%
            {
              fbasis = fda::create.bspline.basis(c(0, I), mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(1:I, x[[j]], fbasis)
              switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*log(I),
                     "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*(log(I) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i],
                     "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          K = numeric(n)
          for(jj in 1:n) K[jj] = mmm[which(v[jj,] == min(v[jj,]))]
          if(commonK == "mode"){
            temp = table(as.vector(K))
            KK[j] = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(commonK == "min"){ KK[j] = min(K) }
          if(commonK == "max"){ KK[j] = max(K) }
          if(commonK == "mean"){ KK[j] = floor(mean(K)) }
        }
      }else{
        if(length(int) != 2){ stop("argument int must be of length two") }
        for(j in 1:p){
          if(parallel.method == "parallel.method0"){
            v = matrix(0, nrow = n, ncol = length(mmm))
            for(i in 1:length(mmm)){
              fbasis = fda::create.bspline.basis(rangeval = c(int[1], int[2]), nbasis = mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(seq(int[1], int[2], length = I), x[[j]], fbasis)
              v[, i] = switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*log(I),
                              "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*(log(I) + 2*gamma.eBIC*log(maxK)),
                              "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i],
                              "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          if(parallel.method == "parallel.method1"){
            v = foreach(i = 1:length(mmm), .combine = cbind) %dopar%
            {
              fbasis = fda::create.bspline.basis(rangeval = c(int[1], int[2]), nbasis = mmm[i], norder = norder)
              fdata = fda::smooth.basisPar(seq(int[1], int[2], length = I), x[[j]], fbasis)
              switch(criterion, "BIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*log(I),
                     "eBIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + mmm[i]*(log(I) + 2*gamma.eBIC*log(maxK)),
                     "AIC" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i],
                     "AICc" = I*log(I*(1-fdata$df/I)^2*fdata$gcv/I) + 2*mmm[i] + (2 * mmm[i] * (mmm[i] + 1))/(n - mmm[i] - 1))
            }
          }
          K = numeric(n)
          for(jj in 1:n) K[jj] = mmm[which(v[jj,] == min(v[jj,]))]
          if(commonK == "mode"){
            temp = table(as.vector(K))
            KK[j] = min(as.numeric(names(temp)[temp == max(temp)]))
          }
          if(commonK == "min"){ KK[j] = min(K) }
          if(commonK == "max"){ KK[j] = max(K) }
          if(commonK == "mean"){ KK[j] = floor(mean(K)) }
        }
      }
    }else{ KK = rep(maxK, p) }
    KM = max(KK)

    data.fd = vector("list", p)
    if(missing(int)){
      fbasisKM = fda::create.bspline.basis(c(0, I), KM, norder = norder)
      for(i in 1:p){
        fbasis = fda::create.bspline.basis(c(0, I), KK[i], norder = norder)
        data.fd[[i]] = t(fda::Data2fd(1:I, x[[i]], fbasis)$coefs)
      }
    }else{
      if(length(int) != 2){ stop("argument int must be of length two") }
      fbasisKM = fda::create.bspline.basis(c(int[1], int[2]), KM, norder = norder)
      for(i in 1:p){
        fbasis = fda::create.bspline.basis(c(int[1], int[2]), KK[i], norder = norder)
        data.fd[[i]] = t(fda::Data2fd(seq(int[1], int[2], length = I), x[[i]], fbasis)$coefs)
      }
    }

    data.ALFA = vector("list", n)
    for(i in 1:l){
      temp = vector("list", n.i[i])
      for(j in 1:n.i[i]){
        temp[[j]] = matrix(rep(0, p*KM), nrow = p)
        for(h in 1:p) temp[[j]][h, 1:KK[h]] =
            data.fd[[h]][group.label == group.label0[i],][j,]
      }
      data.ALFA[group.label == group.label0[i]] = temp
    }

    bspline.cross.prod.matrix = fda::inprod(fbasisKM, fbasisKM)
    AA = 0; BB = 0; CC = 0
    for(i in 1:l){
      z1 = data.ALFA[group.label == group.label0[i]]
      z2 = lapply(z1, t)
      BB1 = 0
      for(j1 in 1:n.i[i]){
        AA = AA + z1[[j1]] %*% bspline.cross.prod.matrix %*% z2[[j1]]
        for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% bspline.cross.prod.matrix %*% z2[[j2]]
        for(m in 1:l){
          z3 = lapply(data.ALFA[group.label == group.label0[m]], t)
          for(j3 in 1:n.i[m]) CC = CC + z1[[j1]] %*% bspline.cross.prod.matrix %*% z3[[j3]]
        }
      }
      BB = BB + BB1/n.i[i]
    }
    CC = CC/n
    Y = solve(AA - CC)
    W.0 = det(AA - BB)
    LH.0 = sum(diag((BB - CC) %*% solve(AA - BB)))
    P.0 = sum(diag((BB - CC) %*% Y))
    R.0 = max(Re(eigen((BB - CC) %*% solve(AA - BB))$values))

    W.P = numeric(B); LH.P = numeric(B); P.P = numeric(B); R.P = numeric(B)
    if(parallel.method == "parallel.method0"){
      for(jP in 1:B){
        BB = 0
        Perm = sample(1:n)
        for(i in 1:l){
          if(i == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(i-1)])+1):sum(n.i[1:i])] }
          z1 = data.ALFA[l1]
          z2 = lapply(z1, t)
          BB1 = 0
          for(j1 in 1:n.i[i]) for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% bspline.cross.prod.matrix %*% z2[[j2]]
          BB = BB + BB1/n.i[i]
        }
        W.P[jP] = det(AA - BB)
        LH.P[jP] = sum(diag((BB - CC) %*% solve(AA - BB)))
        P.P[jP] = sum(diag((BB - CC) %*% Y))
        R.P[jP] = max(Re(eigen((BB - CC) %*% solve(AA - BB))$values))
      }
      W = W.0/det(AA - CC); pvalue.W = mean(W.P <= W.0); LH = LH.0; pvalue.LH = mean(LH.P >= LH.0)
      P = P.0; pvalue.P = mean(P.P >= P.0); R = R.0; pvalue.R = mean(R.P >= R.0)
    }
    if(parallel.method == "parallel.method1"){
      rs.perm = foreach(jP = 1:B, .combine = rbind) %dopar%
      {
        BB = 0
        Perm = sample(1:n)
        for(i in 1:l){
          if(i == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(i-1)])+1):sum(n.i[1:i])] }
          z1 = data.ALFA[l1]
          z2 = lapply(z1, t)
          BB1 = 0
          for(j1 in 1:n.i[i]) for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% bspline.cross.prod.matrix %*% z2[[j2]]
          BB = BB + BB1/n.i[i]
        }
        c(det(AA-BB), sum(diag((BB - CC) %*% solve(AA - BB))),
          sum(diag((BB - CC) %*% Y)), max(Re(eigen((BB - CC) %*% solve(AA - BB))$values)))
      }
      W = W.0/det(AA - CC); pvalue.W = mean(rs.perm[,1] <= W.0); LH = LH.0; pvalue.LH = mean(rs.perm[,2] >= LH.0)
      P = P.0; pvalue.P = mean(rs.perm[,3] >= P.0); R = R.0; pvalue.R = mean(rs.perm[,4] >= R.0)
      parallel::stopCluster(cl)
    }
  }
  if(basis == "own"){
    data.fd = lapply(own.basis, t)
    n = nrow(data.fd[[1]])
    p = length(data.fd)
    if(n != length(group.label)){ stop("the number of observations is not equal to the number of elements in vector of labels") }
    KK = numeric(p)
    for(ii in 1:p) KK[ii] = ncol(data.fd[[ii]])
    KM = max(KK)

    if(KM != nrow(own.cross.prod.mat)){ stop("invalid argument own.cross.prod.mat") }
    if(KM != ncol(own.cross.prod.mat)){ stop("invalid argument own.cross.prod.mat") }

    data.ALFA = vector("list", n)
    for(i in 1:l){
      temp = vector("list", n.i[i])
      for(j in 1:n.i[i]){
        temp[[j]] = matrix(rep(0, p*KM), nrow = p)
        for(h in 1:p) temp[[j]][h, 1:KK[h]] =
            data.fd[[h]][group.label == group.label0[i],][j,]
      }
      data.ALFA[group.label == group.label0[i]] = temp
    }

    AA = 0; BB = 0; CC = 0
    for(i in 1:l){
      z1 = data.ALFA[group.label == group.label0[i]]
      z2 = lapply(z1, t)
      BB1 = 0
      for(j1 in 1:n.i[i]){
        AA = AA + z1[[j1]] %*% own.cross.prod.mat %*% z2[[j1]]
        for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% own.cross.prod.mat %*% z2[[j2]]
        for(m in 1:l){
          z3 = lapply(data.ALFA[group.label == group.label0[m]], t)
          for(j3 in 1:n.i[m]) CC = CC + z1[[j1]] %*% own.cross.prod.mat %*% z3[[j3]]
        }
      }
      BB = BB + BB1/n.i[i]
    }
    CC = CC/n
    Y = solve(AA - CC)
    W.0 = det(AA - BB)
    LH.0 = sum(diag((BB - CC) %*% solve(AA - BB)))
    P.0 = sum(diag((BB - CC) %*% Y))
    R.0 = max(Re(eigen((BB - CC) %*% solve(AA - BB))$values))

    W.P = numeric(B); LH.P = numeric(B); P.P = numeric(B); R.P = numeric(B)
    if(parallel.method == "parallel.method0"){
      for(jP in 1:B){
        BB = 0
        Perm = sample(1:n)
        for(i in 1:l){
          if(i == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(i-1)])+1):sum(n.i[1:i])] }
          z1 = data.ALFA[l1]
          z2 = lapply(z1, t)
          BB1 = 0
          for(j1 in 1:n.i[i]) for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% own.cross.prod.mat %*% z2[[j2]]
          BB = BB + BB1/n.i[i]
        }
        W.P[jP] = det(AA - BB)
        LH.P[jP] = sum(diag((BB - CC) %*% solve(AA - BB)))
        P.P[jP] = sum(diag((BB - CC) %*% Y))
        R.P[jP] = max(Re(eigen((BB - CC) %*% solve(AA - BB))$values))
      }
      W = W.0/det(AA - CC); pvalue.W = mean(W.P <= W.0); LH = LH.0; pvalue.LH = mean(LH.P >= LH.0)
      P = P.0; pvalue.P = mean(P.P >= P.0); R = R.0; pvalue.R = mean(R.P >= R.0)
    }
    if(parallel.method == "parallel.method1"){
      rs.perm = foreach(jP = 1:B, .combine = rbind) %dopar%
      {
        BB = 0
        Perm = sample(1:n)
        for(i in 1:l){
          if(i == 1){ l1 = Perm[1:n.i[1]] }else{ l1 = Perm[(sum(n.i[1:(i-1)])+1):sum(n.i[1:i])] }
          z1 = data.ALFA[l1]
          z2 = lapply(z1, t)
          BB1 = 0
          for(j1 in 1:n.i[i]) for(j2 in 1:n.i[i]) BB1 = BB1 + z1[[j1]] %*% own.cross.prod.mat %*% z2[[j2]]
          BB = BB + BB1/n.i[i]
        }
        c(det(AA-BB), sum(diag((BB - CC) %*% solve(AA - BB))),
          sum(diag((BB - CC) %*% Y)), max(Re(eigen((BB - CC) %*% solve(AA - BB))$values)))
      }
      W = W.0/det(AA - CC); pvalue.W = mean(rs.perm[,1] <= W.0); LH = LH.0; pvalue.LH = mean(rs.perm[,2] >= LH.0)
      P = P.0; pvalue.P = mean(rs.perm[,3] >= P.0); R = R.0; pvalue.R = mean(rs.perm[,4] >= R.0)
      parallel::stopCluster(cl)
    }
  }
  if(basis == "own"){
    result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                  P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                  group.label = group.label, B = B,
                  parallel = parallel, nslaves = nslaves,
                  basis = basis, own.basis = own.basis,
                  own.cross.prod.mat = own.cross.prod.mat)
  }else{
    if(criterion == "NO"){
      if(missing(int)){
        if(basis == "Fourier"){
          result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                        P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                        data = x, group.label = group.label, B = B,
                        parallel = parallel, nslaves = nslaves,
                        basis = basis, criterion = criterion, Km = KK, KM = KM, maxK = maxK)
        }
        if(basis == "b-spline"){
          result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                        P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                        data = x, group.label = group.label, B = B,
                        parallel = parallel, nslaves = nslaves,
                        basis = basis, criterion = criterion, Km = KK, KM = KM, maxK = maxK,
                        norder = norder)
        }
      }else{
        if(basis == "Fourier"){
          result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                        P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                        data = x, group.label = group.label, int = int, B = B,
                        parallel = parallel, nslaves = nslaves,
                        basis = basis, criterion = criterion, Km = KK, KM = KM, maxK = maxK)
        }
        if(basis == "b-spline"){
          result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                        P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                        data = x, group.label = group.label, int = int, B = B,
                        parallel = parallel, nslaves = nslaves,
                        basis = basis, criterion = criterion, Km = KK, KM = KM, maxK = maxK,
                        norder = norder)
        }
      }
    }else{
      if(missing(int)){
        if(criterion == "eBIC"){
          if(basis == "Fourier"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK,
                          gamma.eBIC = gamma.eBIC)
          }
          if(basis == "b-spline"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK, norder = norder,
                          gamma.eBIC = gamma.eBIC)
          }
        }else{
          if(basis == "Fourier"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK)
          }
          if(basis == "b-spline"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK, norder = norder)
          }
        }
      }else{
        if(criterion == "eBIC"){
          if(basis == "Fourier"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, int = int, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK,
                          gamma.eBIC = gamma.eBIC)
          }
          if(basis == "b-spline"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, int = int, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK, norder = norder,
                          gamma.eBIC = gamma.eBIC)
          }
        }else{
          if(basis == "Fourier"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, int = int, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK)
          }
          if(basis == "b-spline"){
            result = list(W = W, pvalueW = pvalue.W, LH = LH, pvalueLH = pvalue.LH,
                          P = P, pvalueP = pvalue.P, R = R, pvalueR = pvalue.R,
                          data = x, group.label = group.label, int = int, B = B,
                          parallel = parallel, nslaves = nslaves,
                          basis = basis, criterion = criterion, commonK = commonK,
                          Km = KK, KM = KM, minK = minK, maxK = maxK, norder = norder)
          }
        }
      }
    }
  }
  class(result) = "fmanovaptbfr"
  return(result)
}
