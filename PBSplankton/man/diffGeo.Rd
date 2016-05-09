\name{diffGeo}
\alias{diffGeo}
\title{
  Difference a Geospatial Object
}
\description{
  Transform a geographic object from a map coordinate system into 
  Cartesian coordinates (units in kilometres) with an origin at or
  near zero or centered on the mean Cartesian coordinate.
}
\usage{
diffGeo(dat, offset.from.zero=0.1, xfld="X", yfld="Y", 
   center=TRUE, replace=FALSE, centroid=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{a geospatial object like a \pkg{PBSmapping} PolySet or EventData set.}
  \item{offset.from.zero}{numeric scalar: essentially an offset from a coordinate origin = 0.}
  \item{xfld}{character string specifying a field in the \code{dat} that contains longitude values.}
  \item{yfld}{character string specifying a field in the \code{dat} that contains latitude values.}
  \item{center}{logical: if \code{TRUE}, centre the data on the mean or a user-specified coordinate.}
  \item{replace}{logical: if \code{TRUE}, replace the original (X,Y) values with the transformed coordinates.}
  \item{centroid}{numeric vector specifying an X- and Y-value on which to centre the data.}
}
\details{
  Input data are returned with additional fields \code{c("Xe", "Yn")},
  arbitrarily named to reflect eastings and northings.
  These fields are the original map coordinates (either \code{"LL"} or
  \code{"UTM"}) transformed to Cartesian coordinates (distance in km)
  from an origin (lower left, at or near zero).

  If the argument \code{center=TRUE}, then additional fields 
  \code{c("cXe", "cYn")} are created, which are the Cartesian 
  coordinates centred on the mean Cartesian coordinate or a user-
  specified map coordinate (internally transformed to 
  Cartesian by the function).

  If \code{replace=TRUE}, then the original map coordinates are copied
  to new fields called \code{"X.orig"} and \code{"Y.orig"}, and the 
  Cartesian coordinated pair -- either \code{c("Xe", "Yn")} or 
  \code{c("cXe", "cYn")}, depending on whether \code{center} is 
  \code{FALSE} or \code{TRUE}, respectively -- is copied to the fields
  \code{c("X","Y")} so that \pkg{PBSmapping} functions can be applied 
  directly.
  
  The projection of the transformed data object is changed from the 
  input object's projection (either \code{"LL"} or \code{"UTM"}) to 
  \code{projection=1}.
}
\value{
  Input data object with additional fields \code{c("Xe", "Yn", "cXe", "cYn")}.
}
\references{
  \href{http://www.movable-type.co.uk/scripts/latlong.html}{Calculate distance, bearing and more between Latitude/Longitude points}
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC
  
  Trevor D. Davies, Stock Assessment Scientist,\cr
  Min. Forests, Lands & Natural Resource Operations, BC Fish & Wildlife Branch, Victoria BC
}
\note{
  This function is used primarily to transform map coordinate data
  for use as an explanator in a GAM (generalized additive model).
}
\seealso{
  \code{\link[PBSmapping]{calcGCdist}},
  \code{\link[PBSplankton]{calcChloro}}
}
\examples{
local(envir=.PBSptonEnv,expr={
  data(nepacLL,surveyData) ## PBSmapping
  coast=nepacLL; survey=surveyData

  ## Need to clip big PolySet to region of the smaller dataset
  ## to avoid radian-to-distance calculation smearing.
  ##----------------------------------------------------------
  survLL  = survey
  coastLL = clipPolys(coast, xlim=range(survLL$X), ylim=range(survLL$Y))

  ## Better to convert the datasets first to UTM coordinates.
  ## Use the same zone when converting using `convUL'.
  ##----------------------------------------------------------
  attr(coast,"zone") = attr(survey, "zone") = 9
  coastUTM = convUL(coast)
  survUTM  = convUL(survey)

  for (i in c("LL","UTM")) {
    sdata = get(paste0("surv",i))
    mdata = get(paste0("coast",i))
    sData = diffGeo(sdata, replace=TRUE)
    mData = diffGeo(mdata, replace=TRUE, 
      centroid=attributes(sData)$centroid$old)
    par(mar=c(4,3,1.2,0.5), mgp=c(1.6,0.5,0))
    plotPoints(sData, projection=1, col="red",pch=20, cex=1.2, plt=NULL)
    abline(v=0, h=0, col="gainsboro")
    addPolys(mData, col=lucent("gainsboro",0.5), border="grey")
  }
})
}
\keyword{IO}
\keyword{algebra}

