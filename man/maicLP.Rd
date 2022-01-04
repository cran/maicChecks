\name{maicLP}
\alias{maicLP}
\title{Checks if AD is within the convex hull of IPD using lp-solve}
\description{Checks if AD is within the convex hull of IPD using lp-solve}
\usage{
maicLP(ipd, ad)
}
\arguments{
  \item{ipd}{a dataframe with n row and p coln, where n is number of subjects and p is the number of variables used in matching.
}
  \item{ad}{a dataframe with 1 row and p coln. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
}
}
\value{
  \item{lp.check }{0 = AD is inside IPD, and MAIC can be conducted; 2 = otherwise}
}
\references{Glimm & Yau (2021). Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons. arXiv 2108.01896.
}
\author{Lillian Yau
}
\examples{
## eAD[1,] is the scenario A in the reference paper,
## i.e. when AD is within IPD convex hull
maicLP(eIPD, eAD[1,2:3])

## eAD[3,] is the scenario C in the reference paper,
## i.e. when AD is outside IPD convex hull
maicLP(eIPD, eAD[3,2:3])
}
