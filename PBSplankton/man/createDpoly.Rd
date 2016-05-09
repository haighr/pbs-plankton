\name{createDpoly}
\alias{createDpoly}
\title{
  Creat a Depth Polygon
}
\description{
  Create a polygon (\pkg{PBSmapping} \code{PolySet}) from one depth contour
  along the BC coast.
}
\usage{
createDpoly(isob=1000, corner=c("TR","BR"), 
   xbox=c(-134,-134,-124.5,-124.5), ybox=c(48,54.4,54.4,48),
   seePlot=FALSE)
}
\arguments{
  \item{isob}{A single isobath (contour) from the set available in
    \pkg{PBSdata}'s \code{\link[PBSdata]{isobath}} object. The options
    only include isobaths from 100m to 1800m at 100m intervals.}
  \item{corner}{Corners needed to close the isobath polygon. Choices:
    \code{"BL"} = bottom left, \code{"TL"} = top left, \code{"TR"} =
    top right, \code{"BR"} = bottom right.}
  \item{xbox}{Four \code{X}-coordinates that outline the bounding box
    of the polygon. Order of vertices follow the order of the corners
    (i.e., clockwise from bottom left).}
  \item{ybox}{Four \code{Y}-coordinates that outline the bounding box
    of the polygon. Order of vertices follow the order of the corners
    (i.e., clockwise from bottom left).}
  \item{seePlot}{logical: if \code{TRUE} draw the coast and the depth
    polygon created by the routine.}
}
\details{
  The diagonal nature of the BC coast means that a depth polygon
  created from an isobath polyline needs corner boundaries.

  Only the longest primary polyline from the isobath set is
  chosen to make the polygon.

  Depending on the complexity of the longest polyline, some 
  experimentation is likely needed before the user chooses the correct
  corners. Even then, additional problems could possibly arise.
  The function has been debugged for the 1000m isobath only.
}
\value{
  A depth polygon is invisibly returned.
  It is also saved as a binary file called \code{isoPoly.rda}.
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC
}
\seealso{
  \code{\link[PBSdata]{isobath}},
  \code{\link[PBSmapping]{clipPolys}},
  \code{\link[PBSmapping]{plotMap}},
  \code{\link[PBSmapping]{addPolys}},
  \code{\link[PBSplankton]{findDpoly}}
}
\examples{
local(envir=.PBSptonEnv,expr={
  depth1000 = createDpoly(seePlot=TRUE)
  print(summary(depth1000))
})
}
\keyword{manip}
\keyword{hplot}

