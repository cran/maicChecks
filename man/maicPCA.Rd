\name{maicPCA}
\alias{maicPCA}
\title{Checks whether AD is outside IPD in PC coordinates
}
\description{Checks whether AD is outside IPD in principal component (PC) coordinates
}
\usage{
maicPCA(ipd, ad)
}
\arguments{
  \item{ipd}{a dataframe with n row and p coln, where n is number of subjects in IPD set and p is the number of variables used in matching.
}
  \item{ad}{a dataframe with 1 row and p coln. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
}
}
\details{When AD is within the IPD PC ranges, AD can still be outside the IPD convex hull in the original scale. On the other hand, if AD is outside the IPD PC ranges, in the original scale AD is for sure outside the IPD convex hull.
}
\value{
Prints a message whether AD is inside or outside IPD PC coordinates. Also returns a ggplot object to be plotted.
  \item{pc.dplot }{dot-plot of AD and IPD both in IPD's PC coordinates}
  \item{pca.check }{0 = AD within the ranges of IPD's PC coordinates; 2 = otherwise}
}
\references{Glimm & Yau (2021). Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons. arXiv 2108.01896.
}
\author{
Lillian Yau
}
\examples{
## eAD[1,] is the scenario A in the reference paper,
## i.e. when AD is perfectly within IPD
a1 <- maicPCA(eIPD, eAD[1,2:3])
a1 ## the dot plots of PC's for IPD and AD

## eAD[3,] is the scenario C in the reference paper,
## i.e. when AD is outside IPD
a3 <- maicPCA(eIPD, eAD[3,2:3])
a3 ## the dot plots of PC's for IPD and AD
}
