## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(TransGraph)
#  # load example data from github repository
#  # Please refer to https://github.com/Ren-Mingyang/example_data_TransGraph for detailed data information
#  load(url("https://github.com/Ren-Mingyang/example_data_TransGraph/raw/main/example.data.GGM.RData"))
#  t.data = example.data.GGM$target.list$t.data
#  t.precision = example.data.GGM$target.list$t.precision
#  A.data = example.data.GGM$A.data
#  A.data.infor = example.data.GGM$A.data.infor
#  
#  # using all auxiliary domains
#  res.trans.weight = trans_precision(t.data, A.data, cov.method="weight")
#  res.trans.opt = trans_precision(t.data, A.data, cov.method="opt")
#  res.trans.size = trans_precision(t.data, A.data, cov.method="size")
#  Theta.trans.weight = res.trans.weight$Theta.hat
#  Theta.trans.opt = res.trans.opt$Theta.hat
#  Theta.trans.size = res.trans.size$Theta.hat
#  Theta.single = res.trans.weight$Theta.hat0  # initial rough estimation via the target domain
#  Theta.single[abs(Theta.single)<0.0001] = 0
#  
#  Evaluation.GGM(Theta.single, t.precision)
#  Evaluation.GGM(Theta.trans.weight, t.precision)
#  Evaluation.GGM(Theta.trans.opt, t.precision)
#  Evaluation.GGM(Theta.trans.size, t.precision)
#  
#  # using informative auxiliary domains
#  res.trans.size.oracle = trans_precision(t.data, A.data.infor, precision.method="CLIME", cov.method="size")
#  Evaluation.GGM(res.trans.size.oracle$Theta.hat, t.precision)
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(TransGraph)
#  library(Tlasso)
#  # load example data from github repository
#  # Please refer to https://github.com/Ren-Mingyang/example_data_TransGraph for detailed data information
#  load(url("https://github.com/Ren-Mingyang/example_data_TransGraph/raw/main/example.data.tensorGGM.RData"))
#  t.data = example.data$t.data
#  A.data = example.data$A.data
#  t.Omega.true.list = example.data$t.Omega.true.list
#  normalize = T
#  
#  K = length(A.data)
#  p.vec = dim(t.data)
#  M = length(p.vec) - 1
#  n = p.vec[M+1]
#  p.vec = p.vec[1:M]
#  tla.lambda = 20*sqrt( p.vec*log(p.vec) / ( n * prod(p.vec) ))
#  A.lambda = list()
#  for (k in 1:K) {
#    A.lambda[[k]] = 20*sqrt( log(p.vec) / ( dim(A.data[[k]])[M+1] * prod(p.vec) ))
#  }
#  
#  # the proposed method
#  res.final = tensor.GGM.trans(t.data, A.data, A.lambda, normalize = normalize)
#  # Tlasso
#  Tlasso.Omega.list = Tlasso.fit(t.data, lambda.vec = tla.lambda, norm.type = 1+as.numeric(normalize))
#  
#  # summary
#  i.Omega = as.data.frame(t(unlist(est.analysis(res.final$Omega.list, t.Omega.true.list))))
#  i.Omega.diff = as.data.frame(t(unlist(est.analysis(res.final$Omega.list.diff, t.Omega.true.list))))
#  i.Tlasso = as.data.frame(t(unlist(est.analysis(Tlasso.Omega.list, t.Omega.true.list))))
#  i.Omega.diff     # proposed.v
#  i.Omega          # proposed
#  i.Tlasso         # Tlasso
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(TransGraph)
#  # load example data from github repository
#  # Please refer to https://github.com/Ren-Mingyang/example_data_TransGraph for detailed data information
#  load(url("https://github.com/Ren-Mingyang/example_data_TransGraph/raw/main/example.data.DAG.RData"))
#  t.data = example.data.DAG$target.DAG.data$X
#  true_adjace = example.data.DAG$target.DAG.data$true_adjace
#  A.data = example.data.DAG$auxiliary.DAG.data$X.list.A
#  
#  # transfer method
#  res.trans = trans.local.DAG(t.data, A.data)
#  # Topological Layer method-based single-task learning (JLMR, 2022)
#  res.single = TLLiNGAM(t.data)
#  
#  Evaluation.DAG(res.trans$B, true_adjace)$Eval_result
#  Evaluation.DAG(res.single$B, true_adjace)$Eval_result
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(TransGraph)
#  # load example data from github repository
#  # Please refer to https://github.com/Ren-Mingyang/example_data_TransGraph for detailed data information
#  load(url("https://github.com/Ren-Mingyang/example_data_TransGraph/raw/main/example.data.singleDAG.RData"))
#  true_adjace = example.data.singleDAG$true_adjace
#  t.data = example.data.singleDAG$X
#  res.single = TLLiNGAM(t.data)
#  Evaluation.DAG(res.single$B, true_adjace)$Eval_result
#  

