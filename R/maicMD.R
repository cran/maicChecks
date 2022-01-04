#
# mahalonobis distance, manually
#
maicMD <- function(ipd, ad, n.ad = Inf) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  md <- y <- NULL
  # ipd and ad must have identical coln's in identical order
  # ipd must be a dataframe
  if (is.data.frame(ipd) == FALSE) {
    ipd <- as.data.frame(ipd)
  }
  # ad must be a row vector
  if (is.data.frame(ad) == FALSE) {
    ad <- data.frame(ad)
    if(dim(ad)[1] > 1) {
      ad <- t(ad)
    }
  }
  # center ipd
  zipd <- as.matrix(sweep(ipd, 2, apply(ipd, 2, mean)))
  # covariance of ipd
  zpz1 <- (dim(ipd)[1] - 1) * solve(t(zipd)%*%(zipd))
  # m-dist. for ipd
  md.ipd <- diag(zipd %*% zpz1 %*% t(zipd))
  # center ad
  zad1 <- as.matrix(ad - apply(ipd, 2, mean))
  # m-dist of ad
  if(n.ad == Inf){
    md.ad <- zad1 %*% zpz1 %*% t(zad1)
  } else{
    md.ad <- n.ad / (dim(ipd)[1] + n.ad) *
      zad1 %*% zpz1 %*% t(zad1)
  }
  if(md.ad > max(md.ipd)) {
    md.check <- 0
  ##  cat('\nM-distance check: AD has the largest Mahalanobis distance to the IPD center.\nTherefore, in original scale AD is *outside* of the IPD convex hull.\n\n')
  } else {
    if(md.ad <= max(md.ipd)){
      md.check <- 2
   ##   cat('\nM-distance check: AD does not have the largest Mahalanobis distance to the IPD center.\n\n')
    }
  }
  #
  # plot
  #
  ipd.plot <- data.frame(md = md.ipd,
                         y = rep(1, length(md.ipd)))
  ad.plot <- data.frame(md = md.ad, y = 1)
  colnames(ad.plot) <- c('md', 'y')
  # coordinate ratio
  # coo.r <- max(md.ipd) / 1.5
  #
  dplot <- ggplot(data = ipd.plot,
                  aes(x = md, y = y),
  ) +
    theme_bw() +
    ylab("") +
    xlab('Mahalanobis distance') +
    geom_point(shape = 1, color = "grey60", size = 2) +
    geom_point(data = ad.plot,
               aes(x = md, y = y),
               shape = 16,
               size = 2.5
    ) +
    theme(panel.grid = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) #+
    #coord_fixed(coo.r)
  #
  return(list(##md.ipd = md.ipd,
              ##md.ad = md.ad,
              md.plot = dplot,
              md.check = md.check)
         )
}

