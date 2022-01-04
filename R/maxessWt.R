maxessWt <- function(ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  ipd.n <- nrow(ipd)
  ones  <- rep(1, ipd.n)
  ipd0  <- data.frame(cbind(ipd, ones))
  ad0   <- matrix(c(ad, 1, rep(0, nrow(ipd0))), nrow = 1)
  p     <- length(ad) # number of variables used in matching
  #
  ## ipd0 serve as the constraint = Amat
  Amat  <- as.matrix(ipd0) # this includes a row of 1's
  ## optimization object: a quadratic function
  Dmat  <- diag(ipd.n)
  Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat)))
  dvec  <- rep(0, ipd.n) # not in the lp-solve
  ## the right hand side
  bvec <- ad0 # include the 1 already
  #
  x1 <- quadprog::solve.QP(Dmat = Dmat,
                           dvec = dvec,
                           Amat = Amat0,
                           bvec = bvec,
                           meq  = p)
  #
  # weights scaled to total number of patients in ipd
  ipd.wts.me <- x1[["solution"]] * ipd.n
  # ess of new weight
  ipd.ess.me <- round(sum(ipd.wts.me)^2 / sum(ipd.wts.me^2), 1)
  ipd.wtsumm.me <- colMeans(ipd * ipd.wts.me)
  ##
  return(list(maxess.wt = ipd.wts.me,
              ipd.ess = ipd.ess.me,
              ipd.wtsumm = ipd.wtsumm.me))
}
