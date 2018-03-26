fmanova.trp = function(x, group.label, k = 30, projection = c("GAUSS", "BM"), 
                       permutation = FALSE, B = 1000, parallel = FALSE, nslaves = NULL) 
{
  n = ncol(x[[1]])
  p = length(x)
  I = nrow(x[[1]])
  if (n != length(group.label)) {
    stop("the number of observations is not equal to the number of elements in vector of labels")
  }
  if (!is.logical(permutation)) {
    stop("argument permutation is not logical")
  }
  if (B < 1) {
    stop("invalid number of permutations (B)")
  }
  if (any(k < 1)) {
    stop("invalid number of projections (k)")
  }
  if (is.null(k)) {
    stop("argument k is not a vector of length greater than or equal to one")
  }
  projection = match.arg(projection)
  if (any(is.na(group.label))) {
    stop("argument group.label can not contain NA values")
  }
  if (!parallel) {
    parallel.method = "parallel.method0"
  }
  else {
    if (!("doParallel" %in% rownames(installed.packages()))) {
      stop("Please install package 'doParallel'")
    }
    requireNamespace("foreach", quietly = TRUE)
    nlp = parallel::detectCores()
    if (is.null(nslaves)) {
      if (nlp >= 2) {
        nslaves = nlp
        parallel.method = "parallel.method1"
      }
      else {
        parallel.method = "parallel.method0"
      }
    }
    else {
      if (nlp >= 2) {
        if (nslaves >= 2) {
          nslaves = nslaves
        }
        else {
          nslaves = nlp
        }
        parallel.method = "parallel.method1"
      }
      else {
        parallel.method = "parallel.method0"
      }
    }
  }
  manova.statistics.quick = function(data, group.label, n, 
                                     p) {
    data = as.matrix(data)
    group.label0 = unique(group.label)
    l = length(group.label0)
    n.i = numeric(l)
    for (i in 1:l) n.i[i] = sum(group.label == group.label0[i])
    gmeans = matrix(colMeans(data), nrow = p, ncol = 1)
    xmeans = matrix(0, nrow = p, ncol = l)
    for (ii in 1:l) xmeans[, ii] = matrix(colMeans(data[group.label == 
                                                          group.label0[ii], ]), nrow = p, ncol = 1)
    H = matrix(0, nrow = p, ncol = p)
    for (iiH in 1:l) H = H + n.i[iiH] * (xmeans[, iiH] - 
                                           gmeans) %*% t(xmeans[, iiH] - gmeans)
    E = matrix(0, nrow = p, ncol = p)
    for (iiE in 1:l) E = E + (n.i[iiE] - 1) * cov(data[group.label == 
                                                         group.label0[iiE], ])
    wart = (eigen(H %*% solve(E)))$values
    return(c(Wilks = prod(1/(1 + wart)), LH = sum(wart), 
             Pillai = sum(wart/(1 + wart)), Roy = wart[1]))
  }
  pvalues = matrix(0, nrow = 4, ncol = length(k))
  all.data.proj = list()
  iik = 0
  modulo = function(z) {
    sqrt(sum(z^2))
  }
  for (ik in k) {
    if (projection == "GAUSS") {
      if (parallel.method == "parallel.method0") {
        ik.data.proj = list()
        W = numeric(ik)
        LH = numeric(ik)
        P = numeric(ik)
        R = numeric(ik)
        for (j in 1:ik) {
          data.proj = matrix(0, nrow = n, ncol = p)
          z = matrix(rnorm(I * p), nrow = p, ncol = I)
          modu = apply(z, 1, modulo)
          z = z/modu
          for (i in 1:p) {
            data.matrix = as.matrix(t(x[[i]]))
            data.proj[, i] = data.matrix %*% z[i, ]
          }
          ik.data.proj[[j]] = data.proj
          model = manova(data.proj ~ as.factor(group.label))
          if (permutation == FALSE) {
            W[j] = summary(model, test = "Wilks")$stats[1, 
                                                        6]
            LH[j] = summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                    6]
            P[j] = summary(model, test = "Pillai")$stats[1, 
                                                         6]
            R[j] = summary(model, test = "Roy")$stats[1, 
                                                      6]
          }
          else {
            Ws = summary(model, test = "Wilks")$stats[1, 
                                                      2]
            LHs = summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                  2]
            Ps = summary(model, test = "Pillai")$stats[1, 
                                                       2]
            Rs = summary(model, test = "Roy")$stats[1, 
                                                    2]
            Wp = numeric(B)
            LHp = numeric(B)
            Pp = numeric(B)
            Rp = numeric(B)
            for (i.perm in 1:B) {
              manovap = manova.statistics.quick(data.proj, 
                                                sample(group.label), n, p)
              Wp[i.perm] = manovap[1]
              LHp[i.perm] = manovap[2]
              Pp[i.perm] = manovap[3]
              Rp[i.perm] = manovap[4]
            }
            W[j] = mean(Re(Wp) <= Ws)
            LH[j] = mean(Re(LHp) >= LHs)
            P[j] = mean(Re(Pp) >= Ps)
            R[j] = mean(Re(Rp) >= Rs)
          }
        }
        iik = iik + 1
        all.data.proj[[iik]] = ik.data.proj
        pvalues[, iik] = c(min(ik * W[order(W)]/1:ik), 
                           min(ik * LH[order(LH)]/1:ik), min(ik * P[order(P)]/1:ik), 
                           min(ik * R[order(R)]/1:ik))
      }
      if (parallel.method == "parallel.method1") {
        cl = parallel::makePSOCKcluster(nslaves)
        doParallel::registerDoParallel(cl)
        rs = foreach(j = 1:ik, .combine = "c") %dopar% 
        {
          data.proj = matrix(0, nrow = n, ncol = p)
          z = matrix(rnorm(I * p), nrow = p, ncol = I)
          modu = apply(z, 1, modulo)
          z = z/modu
          for (i in 1:p) {
            data.matrix = as.matrix(t(x[[i]]))
            data.proj[, i] = data.matrix %*% z[i, ]
          }
          model = manova(data.proj ~ as.factor(group.label))
          if (permutation == FALSE) {
            list(c(summary(model, test = "Wilks")$stats[1, 
                                                        6], summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                                                            6], summary(model, test = "Pillai")$stats[1, 
                                                                                                                                                      6], summary(model, test = "Roy")$stats[1, 
                                                                                                                                                                                             6]), data.proj)
          }
          else {
            Ws = summary(model, test = "Wilks")$stats[1, 
                                                      2]
            LHs = summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                  2]
            Ps = summary(model, test = "Pillai")$stats[1, 
                                                       2]
            Rs = summary(model, test = "Roy")$stats[1, 
                                                    2]
            Wp = numeric(B)
            LHp = numeric(B)
            Pp = numeric(B)
            Rp = numeric(B)
            for (i.perm in 1:B) {
              manovap = manova.statistics.quick(data.proj, 
                                                sample(group.label), n, p)
              Wp[i.perm] = manovap[1]
              LHp[i.perm] = manovap[2]
              Pp[i.perm] = manovap[3]
              Rp[i.perm] = manovap[4]
            }
            list(c(mean(Re(Wp) <= Ws), mean(Re(LHp) >= 
                                              LHs), mean(Re(Pp) >= Ps), mean(Re(Rp) >= 
                                                                               Rs)), data.proj)
          }
        }
        parallel::stopCluster(cl)
        iik = iik + 1
        all.data.proj[[iik]] = rs[seq(2, 2 * ik, by = 2)]
        rs1 = rs[seq(1, 2 * ik - 1, by = 2)]
        rs = matrix(0, nrow = ik, ncol = 4)
        for (jj in 1:ik) rs[jj, ] = rs1[[jj]]
        pvalues[, iik] = c(min(ik * (rs[, 1])[order(rs[, 
                                                       1])]/1:ik), min(ik * (rs[, 2])[order(rs[, 2])]/1:ik), 
                           min(ik * (rs[, 3])[order(rs[, 3])]/1:ik), min(ik * 
                                                                           (rs[, 4])[order(rs[, 4])]/1:ik))
      }
    }
    else {
      if (parallel.method == "parallel.method0") {
        ik.data.proj = list()
        W = numeric(ik)
        LH = numeric(ik)
        P = numeric(ik)
        R = numeric(ik)
        for (j in 1:ik) {
          data.proj = matrix(0, nrow = n, ncol = p)
          for (i in 1:p) {
            data.matrix = as.matrix(t(x[[i]]))
            bm.p = cumsum(rnorm(I, mean = 0, sd = 1))/sqrt(I)
            bm.p = bm.p/modulo(bm.p)
            data.proj[, i] = data.matrix %*% as.matrix(bm.p)
          }
          ik.data.proj[[j]] = data.proj
          model = manova(data.proj ~ as.factor(group.label))
          if (permutation == FALSE) {
            W[j] = summary(model, test = "Wilks")$stats[1, 
                                                        6]
            LH[j] = summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                    6]
            P[j] = summary(model, test = "Pillai")$stats[1, 
                                                         6]
            R[j] = summary(model, test = "Roy")$stats[1, 
                                                      6]
          }
          else {
            Ws = summary(model, test = "Wilks")$stats[1, 
                                                      2]
            LHs = summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                  2]
            Ps = summary(model, test = "Pillai")$stats[1, 
                                                       2]
            Rs = summary(model, test = "Roy")$stats[1, 
                                                    2]
            Wp = numeric(B)
            LHp = numeric(B)
            Pp = numeric(B)
            Rp = numeric(B)
            for (i.perm in 1:B) {
              manovap = manova.statistics.quick(data.proj, 
                                                sample(group.label), n, p)
              Wp[i.perm] = manovap[1]
              LHp[i.perm] = manovap[2]
              Pp[i.perm] = manovap[3]
              Rp[i.perm] = manovap[4]
            }
            W[j] = mean(Re(Wp) <= Ws)
            LH[j] = mean(Re(LHp) >= LHs)
            P[j] = mean(Re(Pp) >= Ps)
            R[j] = mean(Re(Rp) >= Rs)
          }
        }
        iik = iik + 1
        all.data.proj[[iik]] = ik.data.proj
        pvalues[, iik] = c(min(ik * W[order(W)]/1:ik), 
                           min(ik * LH[order(LH)]/1:ik), min(ik * P[order(P)]/1:ik), 
                           min(ik * R[order(R)]/1:ik))
      }
      if (parallel.method == "parallel.method1") {
        cl = parallel::makePSOCKcluster(nslaves)
        doParallel::registerDoParallel(cl)
        rs = foreach(j = 1:ik, .combine = "c") %dopar% 
        {
          data.proj = matrix(0, nrow = n, ncol = p)
          for (i in 1:p) {
            data.matrix = as.matrix(t(x[[i]]))
            bm.p = cumsum(rnorm(I, mean = 0, sd = 1))/sqrt(I)
            bm.p = bm.p/modulo(bm.p)
            data.proj[, i] = data.matrix %*% as.matrix(bm.p)
          }
          model = manova(data.proj ~ as.factor(group.label))
          if (permutation == FALSE) {
            list(c(summary(model, test = "Wilks")$stats[1, 
                                                        6], summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                                                            6], summary(model, test = "Pillai")$stats[1, 
                                                                                                                                                      6], summary(model, test = "Roy")$stats[1, 
                                                                                                                                                                                             6]), data.proj)
          }
          else {
            Ws = summary(model, test = "Wilks")$stats[1, 
                                                      2]
            LHs = summary(model, test = "Hotelling-Lawley")$stats[1, 
                                                                  2]
            Ps = summary(model, test = "Pillai")$stats[1, 
                                                       2]
            Rs = summary(model, test = "Roy")$stats[1, 
                                                    2]
            Wp = numeric(B)
            LHp = numeric(B)
            Pp = numeric(B)
            Rp = numeric(B)
            for (i.perm in 1:B) {
              manovap = manova.statistics.quick(data.proj, 
                                                sample(group.label), n, p)
              Wp[i.perm] = manovap[1]
              LHp[i.perm] = manovap[2]
              Pp[i.perm] = manovap[3]
              Rp[i.perm] = manovap[4]
            }
            list(c(mean(Re(Wp) <= Ws), mean(Re(LHp) >= 
                                              LHs), mean(Re(Pp) >= Ps), mean(Re(Rp) >= 
                                                                               Rs)), data.proj)
          }
        }
        parallel::stopCluster(cl)
        iik = iik + 1
        all.data.proj[[iik]] = rs[seq(2, 2 * ik, by = 2)]
        rs1 = rs[seq(1, 2 * ik - 1, by = 2)]
        rs = matrix(0, nrow = ik, ncol = 4)
        for (jj in 1:ik) rs[jj, ] = rs1[[jj]]
        pvalues[, iik] = c(min(ik * (rs[, 1])[order(rs[, 
                                                       1])]/1:ik), min(ik * (rs[, 2])[order(rs[, 2])]/1:ik), 
                           min(ik * (rs[, 3])[order(rs[, 3])]/1:ik), min(ik * 
                                                                           (rs[, 4])[order(rs[, 4])]/1:ik))
      }
    }
  }
  result = list(pvalues = pvalues, data.projections = all.data.proj, 
                data = x, group.label = group.label, k = k, projection = projection, 
                permutation = permutation, B = B, parallel = parallel, 
                nslaves = nslaves)
  class(result) = "fmanovatrp"
  return(result)
}
