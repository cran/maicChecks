maicLP <- function(ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  ## constrain the sum of the weights to 1
  ones  <- rep(1, nrow(ipd))
  ipd   <- data.frame(cbind(ipd, ones))
  ad    <- data.frame(c(ad, 1))
  p     <- ncol(ad)
  ## the ipd serve as the constraint
  f.con <- as.matrix(t(ipd))
  ## a dummy object to be optimized
  f.obj <- rep(0.5, ncol(f.con))
  ## the right hand side is ad
  f.rhs <- ad
  ## direction of constraint
  f.dir <- rep("=", p)
  ## solve
  lp.check <- lpSolve::lp ("max", f.obj, f.con, f.dir, f.rhs)$status
  ##
  return(list(lp.check = lp.check))
}
