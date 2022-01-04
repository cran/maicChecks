maicPCA <- function (ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  pc.o <- NULL
  ##
  ## standardize ipd and ad
  ##
  ipdm <- apply(ipd, 2, mean)
  ipdsd <- apply(ipd, 2, sd)
  zi1 <- sweep(ipd, 2, ipdm) # scaling done in prcomp()
  w <- (ad - ipdm) / ipdsd # standardize ad w.r.t. ipd
  ##
  ## pca
  ##
  pc <- stats::prcomp(zi1, retx = TRUE, scale. = TRUE)
  pc.ipd <- pc$x
  ## ad in ipd's pca coordinates
  pc.ad <- data.frame(t(pc$rotation) %*% t(w))
  ##
  ## check if ad is within all ipd's pc coordinates
  ##
  x <- data.frame(t(apply(pc.ipd, 2, range)))
  colnames(x) <- c('min', 'max')
  if (all(data.table::between(t(pc.ad), x$min, x$max))) {
    pc.check <- 0 ## pca.check <- 0 ... return asks for pc.check
  ##  cat("PCA check: AD within the ranges of IPD's PC coordinates.\n")
  }
  else {
    pc.check <- 2 ## pca.check <- 2 ... return asks for pc.check
  ##  cat("PCA check: AD outside the ranges of IPD's PC coordinates.\n")
  }
  ##
  ## create plot
  ##
  pc.ipd.long <- tidyr::gather(as.data.frame(pc.ipd), pc, x)
  pc.ipd.long$pc.o <- rep(1:dim(ipd)[2], each = dim(ipd)[1])
  ## ad in pca scale
  pc.ad$pc.o <- 1:nrow(pc.ad)
  colnames(pc.ad) <- c("w", "pc.o")
  ##
  ## plot
  ##
  pc.dplot <-
    ggplot(pc.ipd.long,
           aes(x, factor(pc.o))) +
    geom_point(shape = 1, color = "grey60", size = 2) +
    ## vertical line to go thru 0, i.e. ipd centers in pc coordinates
    geom_vline(xintercept = 0,
               color = "gray25",
               size = .5,
               alpha = .5,
               linetype = 'dashed') +
    geom_point(data = pc.ad,
               mapping = aes(w, factor(pc.o)),
               color = "black",
               size = 2.5,
               shape = 17) +
    theme_bw(base_size = 10) +
    scale_y_discrete(breaks = 1:max(pc.ipd.long$pc.o),
                     labels = paste0("PC", 1:max(pc.ipd.long$pc.o))) +
    ylab("") +
    geom_hline(yintercept = seq(1.5,
                                max(pc.ipd.long$pc.o),
                                by = 1),
               color = "gray",
               size = 0.5,
               alpha = 0.5
    ) +
    guides(size = guide_legend("# of obs.")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor.x = element_blank()) +
    xlab("IPD PC values")
  ##
  return(list(pc.dplot = pc.dplot,
              pc.check = pc.check)
         )
}
