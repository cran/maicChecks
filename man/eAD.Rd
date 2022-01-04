\name{eAD}
\alias{eAD}
\docType{data}
\title{Example AD data}
\description{Example AD data}
\usage{data("eAD")}
\format{
  A data frame with 3 observations on the following 3 variables.
  \describe{
    \item{\code{scen}}{corresponds to scenarios A, B, and C in the reference manuscript (Glimm & Yau (2021)). Scenario A is very close to IPD center (see data(ipd)) and is within the IPD convex hull; scenario B is further away from IPD center but otherwise still inside the IPD convex hull; scenario C is outside IPD convex hull.}
    \item{\code{y1}}{a numeric vector}
    \item{\code{y2}}{a numeric vector}
  }
}
\references{Glimm & Yau (2021). Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons. arXiv 2108.01896.
}
\examples{
print(eAD)
}
\keyword{datasets}
