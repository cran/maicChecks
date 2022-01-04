maicT2Test <- function(ipd, ad, n.ad = Inf){
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  ## n.ad is the number of subjects in the study providing ad.
  ## ... n.ad = NULL implies ad is fixed
  ##
  ipd.bar <- colMeans(ipd)
  n.ipd <- nrow(ipd)
  p.var <- ncol(ipd)
  ## cov.matrix of ipdx
  sig.hat <- (1 / (n.ipd - 1)) * t(ipd) %*%
    (diag(n.ipd) - matrix(1/n.ipd,
                          nrow = n.ipd,
                          ncol = n.ipd)) %*%
    as.matrix(ipd)
  ## hotelling's t-squared
  if(n.ad == Inf){ # when ad is fixed
    T.squared <- n.ipd * as.matrix(ipd.bar - ad) %*%
      solve(sig.hat) %*%
      as.matrix(t(ipd.bar - ad))
  } else {
    T.squared <- ((n.ipd * n.ad) / (n.ipd + n.ad)) *
      as.matrix(ipd.bar - ad) %*%
      solve(sig.hat) %*%
      as.matrix(t(ipd.bar - ad))
  }
  ## hotelling's t-squared adjusting the factor to have an F-distn
  T.squared.f <- T.squared * (n.ipd - p.var) / (p.var * (n.ipd - 1))
  ##
  cat('T.sq.f = ', T.squared.f, '\n')
  p.val <- 1 - pf(T.squared.f, p.var, (n.ipd - p.var))
  cat('p.val = ', p.val, '\n')
}
