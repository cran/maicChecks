maicWt <- function(ipd, ad, max.it = 25) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  objfn <- function(a1, X){
    sum(exp(X %*% a1))
  }
  gradfn <- function(a1, X){
    colSums(sweep(X, 1, exp(X %*% a1), "*"))
  }
  ## to take care the case when only 1 variable is matched
  ipd <- as.data.frame(ipd)
  x <- ncol(ipd)
  ipd.n <- nrow(ipd)
  #
  X.EM.1 <- sweep(as.matrix(ipd), 2, as.matrix(ad), '-') # , give.warning = FALSE
  #
  # Estimate $\alpha_2$ (See Philippo 2016)
  #
  op2 <- stats::optim(par = rep(0, x),
               fn = objfn,
               gr = gradfn,
               X = X.EM.1,
               method = "BFGS",
               control = list(maxit = max.it)
  )
  a2 <- op2$par
  wt <- exp(X.EM.1 %*% a2) # weights for each subject in IPD
  wt.rs <- (wt / sum(wt)) * ipd.n # rescaled weights
  #
  ipd.ess <- round(sum(wt.rs)^2 / sum(wt.rs^2), 1)
  ipd.wtsumm <- colMeans(ipd * wt.rs)
  #
  return(list(optim.out = op2,
              maic.wt = wt,
              maic.wt.rs = wt.rs,
              ipd.ess = ipd.ess,
              ipd.wtsumm = ipd.wtsumm))
}
