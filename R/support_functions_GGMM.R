Init_trans_GGMM = function(t.data, lambda.t, M, A.data, lambda.A.list, M.A.vec,
                           initial.selection="K-means", trace=F ){
  # This function will be improved in the next version
  p = dim(t.data)[2]
  res.target = GGMPF(lambda.t, t.data, M, initial.selection=initial.selection, trace = trace)
  t.Theta_hat.array0 = res.target$opt_Theta_hat
  M0.hat = res.target$K_hat
  t.mean = res.target$opt_Mu_hat
  A.mean = t.mean

  K = length(A.data)
  res.aux.list = list()
  for (k in 1:K) {
    res.aux.list[[k]] = res.target
  }

  #### auxiliary covarianve matrices
  # auxiliary pseudo.cov: Bayesian posterior
  A.cov.soft = list()
  nA.vec.soft = c()
  v = 1
  for (k in 1:K) {
    A.data.k = A.data[[k]]
    for (m in M0.hat) {
      A.cov.soft[[v]] = solve(t.Theta_hat.array0[,,m])
      v = v+1
    }
  }
  A.cov.hard = A.cov.soft
  nA.vec.hard = nA.vec.soft

  res = list(t.Theta_hat.array0=t.Theta_hat.array0, M0.hat=M0.hat, t.mean=t.mean,
             A.cov.soft=A.cov.soft, A.cov.hard=A.cov.hard, A.mean=A.mean,
             nA.vec.soft=nA.vec.soft, nA.vec.hard=nA.vec.hard,
             res.target=res.target, res.aux.list=res.aux.list)
  return(res)

}

f.den.vec = function(data, mu, Theta){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: f.den.vec
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            calculate the density function values at each sample point.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ data: n * p matrix, the design matrix.
  ## @ mu1: p * 1 vector, the mean vector.
  ## @ Omega1: p * p matrix, the precision matrix.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ fdensity: The density function values at each sample point.
  ## -----------------------------------------------------------------------------------------------------------------

  p = length(mu)
  fden = as.numeric( (2*pi)^(-p/2) * (det(Theta))^(1/2) * exp(-1/2*diag(t(t(data) - as.numeric(mu)) %*% Theta %*% (t(data) - as.numeric(mu)))) )
  return(fden)
}


