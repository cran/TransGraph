trans.estimate = function(t.data.tran, A.data, A.lambda, A.orac = NULL,
                          t.lambda.int=NULL, adjust.BIC=FALSE, mode.set = NULL,
                          init.method="sepa", init.method.aux="sepa",
                          init.iter=3, normalize = TRUE,
                          theta.algm="cd", cov.select="tensor.prod",
                          c.lam.sepa=20, c.lam.Tlasso=20, cn.lam2=1,
                          inti.the=TRUE){
  
  p.vec = dim(t.data.tran)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]
  if(is.null(mode.set)){
    mode.set = 1:M
  } else {
    mode.set = mode.set
  }
  
  K = length(A.data)
  nA.vec = rep(0, K)
  for (k in 1:K) {
    p.vec.A = dim(A.data[[k]])
    nA.vec[k] = p.vec.A[length(p.vec.A)]
  }
  
  # Initialization
  if(is.null(t.lambda.int)){
    if(init.method == "sepa"){
      t.lambda.tran = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(init.method == "Tlasso"){
      t.lambda.tran = c.lam.Tlasso*sqrt( log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(M==2){
      t.lambda.tran = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
  }else{
    t.lambda.tran = t.lambda.int
  }
  
  init.res = Initial(t.data.tran, t.lambda.tran, A.data, A.lambda,
                     A.orac = A.orac, method = init.method, method.aux = init.method.aux,
                     TT=init.iter, normalize = normalize, mode.set = mode.set)
  init.time = init.res$time
  t.Omega.hat.list = init.res$t.Omega.hat.list
  
  
  if(cov.select=="inverse"){
    # 0 Covariance matrix: directly inverting precision matrix
    S.hat.A.list = init.res$S.hat.A.M0
    S.hat.A.diff.list = init.res$S.hat.A.M.diff0
  }
  if(cov.select=="tensor.prod"){
    # 1 Covariance matrix: multiplication by tensors and precision matrices
    S.hat.A.list = init.res$S.hat.A.M1
    S.hat.A.diff.list = init.res$S.hat.A.M.diff1
  }
  
  t1 = proc.time()
  # using the weight (in Sigma.A.m) determined by the differences
  Theta.hat.diff.list = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    S.hat.A = S.hat.A.diff.list[[mi]]
    Omega.hat0 = t.Omega.hat.list[[m]]
    lam1 = 2*max(apply(abs(Omega.hat0), 2, sum))*sqrt(log(p.vec[m]) / n.da)
    delta.hat = delta.est(S.hat.A, Omega.hat0, lam1=lam1)
    
    h.hat = max(apply(delta.hat, 2, function(x) sum(abs(x))))
    lam2.del = min(h.hat*sqrt(p.vec[m]*log(p.vec[m]) / n.da / prod(p.vec)), h.hat^2)
    lam2.N = sqrt(p.vec[m]*log(p.vec[m]) / sum(nA.vec) / prod(p.vec))
    lambda2 = cn.lam2*min(max(lam2.del, lam2.N), n.da*lam2.N)
    
    Omega.hat00 = Omega.hat0 * inti.the + 0 * (!inti.the)
    Theta.tuning.res = Theta.tuning(lambda2, S.hat.A, delta.hat, Omega.hat00,
                                    n.A=sum(nA.vec), theta.algm=theta.algm, adjust.BIC=adjust.BIC)
    Theta.hat.m = Theta.tuning.res$Theta.hat.m
    Theta.hat.diff.list[[mi]] = Theta.hat.m
  }
  
  ## using the weight (in Sigma.A.m) determined by the sample sizes
  Theta.hat.list = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    S.hat.A = S.hat.A.list[[mi]]
    Omega.hat0 = t.Omega.hat.list[[m]]
    lam1 = 2*max(apply(abs(Omega.hat0), 2, sum))*sqrt(log(p.vec[m]) / n.da)
    delta.hat = delta.est(S.hat.A, Omega.hat0, lam1=lam1)
    
    h.hat = max(apply(delta.hat, 2, function(x) sum(abs(x))))
    lam2.del = min(h.hat*sqrt(p.vec[m]*log(p.vec[m]) / n.da / prod(p.vec)), h.hat^2)
    lam2.N = sqrt(p.vec[m]*log(p.vec[m]) / sum(nA.vec) / prod(p.vec))
    lambda2 = cn.lam2*min(max(lam2.del, lam2.N), n.da*lam2.N)
    
    Omega.hat00 = Omega.hat0 * inti.the + 0 * (!inti.the)
    Theta.tuning.res = Theta.tuning(lambda2, S.hat.A, delta.hat, Omega.hat00,
                                    n.A=sum(nA.vec), theta.algm=theta.algm, adjust.BIC=adjust.BIC)
    Theta.hat.m = Theta.tuning.res$Theta.hat.m
    Theta.hat.list[[mi]] = Theta.hat.m
  }
  t.theta = proc.time() - t1
  
  if(sum(A.orac) > 0){
    ########### oracle auxiliary domains
    if(cov.select=="inverse"){
      S.hat.A.list.o = init.res$S.hat.A.M0.o   # 0 Covariance matrix: directly inverting precision matrix
    }
    if(cov.select=="tensor.prod"){
      S.hat.A.list.o = init.res$S.hat.A.M1.o   # 1 Covariance matrix: multiplication by tensors and precision matrices
    }
    
    Theta.hat.list.o = list()
    for (m in mode.set) {
      mi = match(m, mode.set)
      S.hat.A = S.hat.A.list.o[[mi]]
      Omega.hat0 = t.Omega.hat.list[[m]]
      lam1 = 2*max(apply(abs(Omega.hat0), 2, sum))*sqrt(log(p.vec[m]) / n.da)
      delta.hat = delta.est(S.hat.A, Omega.hat0, lam1=lam1)
      
      lambda2 = cn.lam2*sqrt(p.vec[m]*log(p.vec[m]) / sum(nA.vec[A.orac]) / prod(p.vec))
      Theta.tuning.res = Theta.tuning(lambda2, S.hat.A, delta.hat, Omega.hat0,
                                      n.A=sum(nA.vec[A.orac]), theta.algm=theta.algm, adjust.BIC=adjust.BIC)
      Theta.hat.m = Theta.tuning.res$Theta.hat.m
      Theta.hat.list.o[[mi]] = Theta.hat.m
    }
    res.trans = list(Theta.hat.list=Theta.hat.list, Theta.hat.diff.list=Theta.hat.diff.list,
                     S.hat.A.list=S.hat.A.list, S.hat.A.diff.list=S.hat.A.diff.list,
                     t.Omega.hat.list=t.Omega.hat.list, init.res=init.res,
                     init.time=init.time, theta.time=t.theta,
                     Theta.hat.list.o=Theta.hat.list.o, S.hat.A.list.o=S.hat.A.list.o)
  } else {
    res.trans = list(Theta.hat.list=Theta.hat.list, Theta.hat.diff.list=Theta.hat.diff.list,
                     S.hat.A.list=S.hat.A.list, S.hat.A.diff.list=S.hat.A.diff.list,
                     t.Omega.hat.list=t.Omega.hat.list, init.res=init.res,
                     init.time=init.time, theta.time=t.theta)
  }
  
  return(res.trans)
  
}

Initial = function(t.data, t.lambda, A.data, A.lambda,
                   A.orac = NULL, method = "sepa", method.aux = "sepa",
                   TT=2, normalize = TRUE, mode.set = NULL){
  # Initial: the function calculating initial precision matrices of the target domain
  #          and covariance matrices of the auxiliary domain,
  #          via two alternative methods:
  #          "Tlasso" (PAMI, 2020) & "sepa" (JCGS, 2022)
  
  p.vec = dim(t.data)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]
  K = length(A.data)
  nA.vec = rep(0, K)
  for (k in 1:K) {
    p.vec.A = dim(A.data[[k]])
    nA.vec[k] = p.vec.A[length(p.vec.A)]
  }
  if(is.null(mode.set)){
    mode.set = 1:M
  } else {
    mode.set = mode.set
  }
  
  t0 <- proc.time()
  if(method == "sepa"){
    # Initialization in target domain
    t.Omega.hat.list = Separate.fit(t.data, lambda.vec=t.lambda, normalize = normalize)$Omegahat
  }
  
  if(method == "Tlasso"){
    # Initialization in target domain
    t.Omega.hat.list = Tlasso.fit(t.data, T=TT, lambda.vec = t.lambda, norm.type = 1+as.numeric(normalize))
  }
  
  if(method.aux == "sepa"){
    # Initialization in auxiliary domains
    A.Omega.hat.list = list()
    for (k in 1:K) {
      A.Omega.hat.list[[k]] = Separate.fit(A.data[[k]], lambda.vec=A.lambda[[k]], normalize = normalize)$Omegahat
    }
  }
  if(method.aux == "Tlasso"){
    # Initialization in auxiliary domains
    A.Omega.hat.list = list()
    for (k in 1:K) {
      A.Omega.hat.list[[k]] = Tlasso.fit(A.data[[k]], lambda.vec = A.lambda[[k]], norm.type = 1+as.numeric(normalize))
    }
  }
  
  A.S.hat.list0 = list()
  A.S.hat.list1 = list()
  for (k in 1:K) {
    A.S.hat.list = S.est(A.data[[k]], A.Omega.hat.list[[k]])
    A.S.hat.list0[[k]] = A.S.hat.list$sig0
    A.S.hat.list1[[k]] = A.S.hat.list$sig1
  }
  
  S.hat.A.M = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    S.hat.A.M[[mi]] = diag(p.vec[m]) - diag(p.vec[m])
  }
  # weight determined by the differences
  S.hat.A.M.diff0 = S.hat.A.M
  S.hat.A.M.diff1 = S.hat.A.M
  weight.KM0 = matrix(0, ncol = K, nrow = length(mode.set))
  weight.KM1 = matrix(0, ncol = K, nrow = length(mode.set))
  for (k in 1:K) {
    for (m in mode.set) {
      mi = match(m, mode.set)
      weight.KM0[mi,k] = 1/sum((A.S.hat.list0[[k]][[m]] %*% t.Omega.hat.list[[m]] - diag(p.vec[m]))^2)
      weight.KM1[mi,k] = 1/sum((A.S.hat.list1[[k]][[m]] %*% t.Omega.hat.list[[m]] - diag(p.vec[m]))^2)
    }
  }
  
  wed0 = t(t(weight.KM0)*nA.vec)
  alpha.k.diff0 = wed0 / apply(wed0, 1, sum)
  wed1 = t(t(weight.KM1)*nA.vec)
  alpha.k.diff1 = wed1 / apply(wed1, 1, sum)
  for (k in 1:K) {
    for (m in mode.set) {
      mi = match(m, mode.set)
      S.hat.A.M.diff0[[mi]] = S.hat.A.M.diff0[[mi]] + A.S.hat.list0[[k]][[m]] * alpha.k.diff0[mi,k]
      S.hat.A.M.diff1[[mi]] = S.hat.A.M.diff1[[mi]] + A.S.hat.list1[[k]][[m]] * alpha.k.diff1[mi,k]
    }
  }
  
  
  # weight determined by the sample sizes
  S.hat.A.M0 = S.hat.A.M
  S.hat.A.M1 = S.hat.A.M
  alpha.k = nA.vec/sum(nA.vec)
  for (k in 1:K) {
    for (m in mode.set) {
      mi = match(m, mode.set)
      S.hat.A.M0[[mi]] = S.hat.A.M0[[mi]] + A.S.hat.list0[[k]][[m]] * alpha.k[k]
      S.hat.A.M1[[mi]] = S.hat.A.M1[[mi]] + A.S.hat.list1[[k]][[m]] * alpha.k[k]
    }
  }
  t00 = proc.time() - t0
  
  
  if(sum(A.orac) > 0){
    S.hat.A.M0.o = S.hat.A.M
    S.hat.A.M1.o = S.hat.A.M
    alpha.k.o = nA.vec[A.orac]/sum(nA.vec[A.orac])
    A.S.hat.list0.o = list()
    A.S.hat.list1.o = list()
    for (k in A.orac) {
      ki = match(k, A.orac)
      A.S.hat.list.o = S.est(A.data[[k]], A.Omega.hat.list[[k]])
      A.S.hat.list0.o[[ki]] = A.S.hat.list.o$sig0
      A.S.hat.list1.o[[ki]] = A.S.hat.list.o$sig1
      for (m in mode.set) {
        mi = match(m, mode.set)
        S.hat.A.M0.o[[mi]] = S.hat.A.M0.o[[mi]] + A.S.hat.list0.o[[ki]][[m]] * alpha.k.o[ki]
        S.hat.A.M1.o[[mi]] = S.hat.A.M1.o[[mi]] + A.S.hat.list1.o[[ki]][[m]] * alpha.k.o[ki]
      }
    }
    Init.res = list(t.Omega.hat.list=t.Omega.hat.list, A.Omega.hat.list=A.Omega.hat.list,
                    S.hat.A.M0=S.hat.A.M0, S.hat.A.M1=S.hat.A.M1,
                    S.hat.A.M.diff0=S.hat.A.M.diff0, S.hat.A.M.diff1=S.hat.A.M.diff1,
                    A.S.hat.list0=A.S.hat.list0, A.S.hat.list1=A.S.hat.list1,
                    S.hat.A.M0.o=S.hat.A.M0.o, S.hat.A.M1.o=S.hat.A.M1.o,
                    A.S.hat.list0.o=A.S.hat.list0.o, A.S.hat.list1.o=A.S.hat.list1.o,
                    time = t00)
    return(Init.res)
  } else {
    Init.res = list(t.Omega.hat.list=t.Omega.hat.list, A.Omega.hat.list=A.Omega.hat.list,
                    S.hat.A.M0=S.hat.A.M0, S.hat.A.M1=S.hat.A.M1,
                    S.hat.A.M.diff0=S.hat.A.M.diff0, S.hat.A.M.diff1=S.hat.A.M.diff1,
                    A.S.hat.list0=A.S.hat.list0, A.S.hat.list1=A.S.hat.list1,
                    time = t00)
    return(Init.res)
  }
  
}


Initial.aggr = function(t.data, t.lambda.int=NULL, method = "sepa",
                        cov.select= "tensor.prod", TT=2,
                        c.lam.sepa=20, c.lam.Tlasso=20, normalize = TRUE){
  # Initial.aggr: the function calculating initial covariance matrices of
  #               the target domain for the aggregation step,
  #               via two alternative methods:
  #               "Tlasso" (PAMI, 2020) & "sepa" (JCGS, 2022)
  
  p.vec = dim(t.data)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]
  
  if(is.null(t.lambda.int)){
    if(method == "sepa"){
      t.lambda = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(method == "Tlasso"){
      t.lambda = c.lam.Tlasso*sqrt( log(p.vec) / ( n.da * prod(p.vec) ))
    }
    if(M==2){
      t.lambda = c.lam.sepa*sqrt( p.vec*log(p.vec) / ( n.da * prod(p.vec) ))
    }
  }else{
    t.lambda = t.lambda.int
  }
  
  
  if(method == "sepa"){
    # Initialization in target domain
    t.Omega.hat.list = Separate.fit(t.data, lambda.vec=t.lambda, normalize=normalize)$Omegahat
  }
  
  if(method == "Tlasso"){
    # Initialization in target domain
    t.Omega.hat.list = Tlasso.fit(t.data, T=TT, t.lambda, norm.type = 1+as.numeric(normalize))
  }
  
  t.S.hat.list = S.est(t.data, t.Omega.hat.list)
  
  t.S.hat.list0 = t.S.hat.list$sig0  # 0 Covariance matrix: directly inverting precision matrix
  t.S.hat.list1 = t.S.hat.list$sig1  # 1 Covariance matrix: multiplication by tensors and precision matrices
  
  if(cov.select=="inverse"){
    t.S.hat.list = t.S.hat.list0
  }
  if(cov.select=="tensor.prod"){
    t.S.hat.list = t.S.hat.list1
  }
  
  Init.res = list(t.Omega.hat.list=t.Omega.hat.list,
                  t.S.hat.list=t.S.hat.list)
  
  
  
}

S.est = function(data, Omega.hat.list){
  # S.est: the function calculating the covariance matrix of each mode
  
  p.vec = dim(data)
  M = length(p.vec) - 1
  n.da = p.vec[M+1]
  p.vec = p.vec[-(M+1)]
  
  # 0 Covariance matrix: directly inverting precision matrix
  Omega.hat.list.sqrt = list()
  S.hat.list0 = list()
  for (m in 1:M) {
    Omega.hat.list.sqrt[[m]] = expm::sqrtm(Omega.hat.list[[m]])
    S.hat.list0[[m]] = solve(Omega.hat.list[[m]])
  }
  
  # 1 Covariance matrix: multiplication by tensors and precision matrices
  S.hat.list1 = list()
  for(m in 1:M){
    S.array = array(0,c(p.vec[m],p.vec[m],n.da))
    Omega.hat.list.sqrt.m = Omega.hat.list.sqrt
    Omega.hat.list.sqrt.m[[m]] = diag(p.vec[m])
    for(i in 1:n.da){
      d=0
      eval(parse(text=paste('d=data[',paste(rep(',',M),collapse=''),'i]')))
      Vi = k_unfold( as.tensor(ttl( as.tensor(d) , Omega.hat.list.sqrt.m , ms=1:M)@data) ,m=m)@data
      S.array[,,i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array,c(1,2),mean) * p.vec[m] / prod(p.vec)
    S.hat.list1[[m]] = S.mat
  }
  S.hat.list = list(sig0=S.hat.list0, sig1=S.hat.list1)
  return(S.hat.list)
}

select.1 = function(t.sigma.tilde.list, t.Omega.hat.list, Theta.hat.list, mode.set, symmetric=TRUE){
  p.vec = NULL
  for(j in 1:length(t.Omega.hat.list)){
    p.vec[j] = dim(t.Omega.hat.list[[j]])[1]
  }
  
  Omega.hat.final.list = list()
  Omega.hat.final.sym.list = list()
  W.list = list()
  for (m in mode.set) {
    mi = match(m, mode.set)
    pm = p.vec[m]
    t.sigma.tilde.m = t.sigma.tilde.list[[m]]
    Theta.trans.m = Theta.hat.list[[mi]]
    t.Omega.hat.m = t.Omega.hat.list[[m]]
    
    resid.Omega0 = apply((t.sigma.tilde.m %*% t.Omega.hat.m - diag(pm))^2, 2, sum)
    resid.trans = apply((t.sigma.tilde.m %*% Theta.trans.m - diag(pm))^2, 2, sum)
    w.m = apply(rbind(resid.Omega0, resid.trans),2,which.min)
    Omega.final.m = t( t(t.Omega.hat.m) * c(w.m == 1) + t(Theta.trans.m) * c(w.m == 2) )
    
    Omega.hat.final.list[[mi]] = Omega.final.m
    W.list[[mi]] = w.m # 1: only target; 2: transfer
    
    if(symmetric){
      Omega.hat.final.sym.list[[mi]] = symmetric.mat(Omega.final.m)}
  }
  selec = list(Omega.hat.final.list = Omega.hat.final.list,
               Omega.hat.final.sym.list = Omega.hat.final.sym.list,
               W.list = W.list)
  return(selec)
  
}

Separate.fit = function(t.data, lambda.vec, normalize){
  warning("The sepa method has not been included in the current version of this R package to circumvent code ownership issues.")
}
