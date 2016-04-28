\name{calcChloro}
\alias{calcChloro}
\title{
  Calculate Chlorophyll Blooms
}
\description{
  Smooth satellite chlorophyll data, calculate bloom periods, and 
  estimate total chlorophyll in blooms by integrating under the curve.
}
\usage{
calcChloro(dat, wd=getwd(), xfld="period", yfld="Chlo", xint=0.1,
   span=0.75, diff.required=FALSE, slice=NULL, days.per=8, 
   salmEnd="Jun-15", show=c(fit=TRUE, der1=TRUE, der2=FALSE), 
   slcs.per.page=1, plts.hide=FALSE, pdf=FALSE, outfile="ChloroFits")
}
\arguments{
  \item{dat}{satellite data file (or name of the file) containing
    fields like \code{year}, \code{X} (longitude), \code{Y} (latitude),
    \code{julian} (days 1 to 365 within a year), \code{chl} 
    (cholorophyll \emph{a} in mg/m\eqn{^3}), and \code{region}.}
  \item{wd}{working directory (defaults to user's working directory).}
  \item{xfld}{field that describes time -- e.g., \code{"julian"}.}
  \item{yfld}{field that describes chlorophyll (or some other oceanic measurement).}
  \item{xint}{time interval for the smoothed line (x-coordinates).}
  \item{span}{the \code{loess} \eqn{\alpha}-parameter that controls the degree of smoothing.}
  \item{diff.required}{logical: if \code{TRUE} then Y-data is assumed
    to be cumulative and requires differencing to get raw data.}
  \item{slice}{list object of discrete logical units -- usually years and regions.}
  \item{days.per}{number of days per period, given that satellite data 
    are often summarised over time periods (e.g., weekly, 8-days, etc.).}
  \item{salmEnd}{date when salmon have left the building (LG to clarify).}
  \item{show}{which fits to display -- \code{fit=T}: fit to smoothed data and integration,
    \code{der1=T}: first derivative fit, \code{der2=T}: second derivative fit.}
  \item{slcs.per.page}{number of slices (e.g., year-region) to display per page.}
  \item{plts.hide}{logical: if \code{TRUE}, no plots are displayed but
    the results table is still created.}
  \item{pdf}{logical: if \code{TRUE}, send the plots to a PDF device (file).}
  \item{outfile}{name for the output image PDF and the results table CSV files.}
}
\details{
  This function attempts to identify phytoplankton bloom periods from
  satellite data, which will determine the zooplpankton response and 
  subsequently, the effect on pelagic fish (particularly salmon) 
  production at sea. 
}
\value{
  A list object where each element contains model results for a slice.\cr
  The model results comprise:\cr
  \code{resChlo} -- output from the function \code{\link[PBSplankton]{calcDeriv}};\cr
  \code{predVal} -- predicted values from smoothed function;\cr
  \code{deriVal} -- results from the derivative function;\cr
  \code{xpeak} -- x-coordinates (time) associated with the peaks in chlorophyll;\cr
  \code{ypeak} -- y-coordinates describing the peak values in chlorophyll.
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC.
}
\note{
  Satellite data from NOAA can be processed using the package \pkg{PBSsatellite}.
}
\seealso{
  \code{\link[PBSplankton]{calcDeriv}},
  \code{\link[numDeriv]{grad}},
  \code{\link[stats]{integrate}}
}
\keyword{hplot}
\keyword{math}

