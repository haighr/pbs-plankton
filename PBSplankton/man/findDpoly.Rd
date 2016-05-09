\name{findDpoly}
\alias{findDpoly}

\title{
  Find Events in Depth Polygon
}
\description{
  Find satellite chorophyll events (gridded) in a polygon bounded 
  seaward by a depth contour.
}
\usage{
findDpoly(edata, pdata, zfld="Chl", latBands, 
   seePlot=TRUE, xlim, ylim, onam="filtered")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edata}{Event data -- usually gridded satellite cholorophyll readings.}
  \item{pdata}{Polygon data -- usually bounded seaward by a depth contour.}
  \item{zfld}{Field name in event data that contains chlorophyll readings.}
  \item{latBands}{List object of named latitude limits within \code{pdata},
    should the user wish to identify events within specific latitudinal bands.}
  \item{seePlot}{logical: if \code{TRUE}, plot the results on a map.}
  \item{xlim}{Longitude limits of the map.}
  \item{ylim}{Latitude limits of the map.}
  \item{onam}{Output name (without extension) of a binary results file (\code{*.rda}).}
}
\details{
  Identifies chlorophyll events that fall into the specified polygon.
  Satellite chlorophyll readings are identified as valid if they lie
  in the interval \code{0<=[Chl]<=30}. Missing values are usually 
  coded as \code{-99}.
}
\value{
  A list object of events with at least one valid chlorophyll reading
  that occur within the specified polygon \code{pdata} if no latitude
  bands are specified or within latitude bands if these are specified
  (e.g., \code{latBands = list(hecate=c(52,54), qcsound=c(51,52))}).
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC
}
\seealso{
  \code{\link[PBSplankton]{createDpoly}},
  \code{\link[PBSplankton]{transData}}
}
\examples{
local(envir=.PBSptonEnv,expr={
  ## If you have files `edata' and `pdata' just load them
  exdir = system.file("examples",package="PBSplankton")
  edata = transData(paste0(exdir,"/wee_SEAWIFS.csv"))
  pdata = createDpoly(seePlot=FALSE)
  findDpoly(edata, pdata, seePlot=TRUE, onam="wee_SEAWIFS_filtered",
    latBands=list(hecate=c(52,54), qcsound=c(51,52)),
    xlim=c(-133.5,-128), ylim=c(50.5,54.5))
})
}
\keyword{manip}
\keyword{hplot}

