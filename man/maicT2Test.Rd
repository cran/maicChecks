\name{maicT2Test}
\alias{maicT2Test}
\title{Hotelling's T-square test to check whether maic is needed
}
\description{Conducts Hotelling's T-square test
}
\usage{
maicT2Test(ipd, ad, n.ad = Inf)
}
\arguments{
  \item{ipd}{a dataframe with n row and p coln, where n is number of subjects and p is the number of variables used in matching.
}
  \item{ad}{a dataframe with 1 row and p coln. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
}
  \item{n.ad}{default is Inf assuming \code{ad} is a fixed (known) quantity with infinit accuracy. In most MAIC applications \code{ad} is the sample statistics and \code{n.ad} is known.
}
}
\details{When \code{n.ad} is not Inf, the covariance matrix is adjusted by the factor n.ad/(n.ipd + n.ad)), where n.ipd is nrow(ipd), the sample size of \code{ipd}.
}
\value{
\item{T.sq.f }{the value of the T^2 test statistic}
\item{p.val }{the p-value corresponding to the test statistic. When the p-value is small, matching is necessary.}
}
\references{Glimm & Yau (2021). Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons. arXiv 2108.01896.
}
\author{
Lillian Yau & Ekkehard Glimm
}
\examples{
## eAD[1,] is the scenario A in the reference paper,
## i.e. when AD is perfectly within IPD
maicT2Test(eIPD, eAD[1,2:3])
}
