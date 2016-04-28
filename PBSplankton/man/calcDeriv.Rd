\name{calcDeriv}
\alias{calcDeriv}
\title{
  Calculate Derivatives of Smoothed Curve
}
\description{
  Calculate the first and second derivatives from a loess-smoothed
  curve and integrate the area beneath it.
}
\usage{
calcDeriv(x, y, span = 0.5)
}
\arguments{
  \item{x}{the x-coordinate data (usually time)}
  \item{y}{the y-coordinate data (e.g., cholorophyll)}
  \item{span}{the \code{loess} \eqn{\alpha}-parameter that controls the degree of smoothing.}
}
\details{
  The function smoothes the (x,y) data using the \code{loess} function,
  calculates the first and second derivatives, then integrates under the 
  loess-curve from \code{min(x)} to \code{max(x)}.
}
\value{
  A list object containing:\cr
  \code{fitline} -- the loess fit,\cr
  \code{predFunc} -- the prediction function,\cr
  \code{deriFunc} -- the function to calculate the first derivative,\cr
  \code{hessFunc} -- the function to calculate the second derivative, and \cr
  \code{integratey} -- the integration of y under the loess-smoothed line.
}
\author{
  Rowan Haigh and Wayne Hajas, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC.
}
\seealso{
See \code{\link[numDeriv]{grad}} in the package \pkg{numDeriv} by P. Gilbert and R. Varadhan.\cr
The function \code{calcDeriv} is used by \code{\link[PBSplankton]{calcChloro}}.
}
\keyword{smooth}
\keyword{math}

