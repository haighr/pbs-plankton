\name{PBSplankton-package}
\alias{PBSplankton-package}
\alias{PBSplankton}
\docType{package}
\title{
  Plankton from Samples to Satellites
}
\description{
  Provides tools for exploring and processing plankton data from samples to satellites.
}
\details{
  \href{https://github.com/pbs-software/pbs-plankton}{PBSplankton} provides an
  R interface for algorithms used in biological (and fisheries) oceanography.
  The scope of this package is by no means comprehensive. Many of the functions
  provide a quick way to visualize data, and in some cases perform preliminary
  analyses. Though oriented to users at the Pacific Biological Station (PBS),
  these functions may prove to be useful to users at other locales. The User's
  Guide (\emph{currently incomplete}) is organised into sections that loosely classify
  the functions by theme -- (1) Utility, (2) Temporal, and (3) Satellite. 
  Within each section, the functions are described alphabetically.

  Additional functions geared toward fisheries can be found in the package
  \href{https://github.com/pbs-software/pbs-tools}{PBStools}, upon which 
  \code{PBSplankton} depends. In turn, \code{PBStools} depends heavily on two other
  R packages: \href{https://github.com/pbs-software/pbs-mapping}{PBSmapping}
  and \href{https://github.com/pbs-software/pbs-modelling}{PBSmodelling}.
  \code{PBSplankton} has a dedicated temporary working environment called 
  \code{.PBSptonEnv}. Accessor functions called \code{ptget}, \code{ptcall},
  and \code{ptprint} enable the user to get, call, and print objects from the
  \code{.PBSptonEnv} environment, while \code{ptput} puts (writes) objects to
  \code{.PBSptonEnv}.
}
\section{PBSsatellite}{
  \href{http://roam.macewan.ca/islandora/object/gm\%3A1451}{PBSsatellite} is a
  collaborative effort between the Pacific Biological Station (PBS) and MacEwan
  University to design software that simplifies the extraction and statistical
  analysis of gridded satellite data. This software uses \code{PBSmapping}, an
  existing R software package, to produce plots. Additionally, users are able
  to produce their own functionality and data interpretation algorithms and
  incorporate these functions into the package's analysis methods. This provides
  users with a platform that simplifies the complex problem of satellite data
  analysis and visualization while providing the ability for personal research
  and method integration.
  
  \emph{Note:} At this stage, \code{PBSplankton} does not take full advantage
  of the functionality in \code{PBSsatellite}.
}
\section{Integrating Chlorophyll}{

\bold{Step 1.} -- Transpose tabulated geo-referenced data from a file into a
  data frame. The routine attempts to automatically detect the longitude and
  latitude fields.
  See \code{\link[PBSplankton]{transData}} for details.
\preformatted{
## If you have the annual data binary file, just load it;
## otherwise load from an ascii file:
exdir = system.file("examples",package="PBSplankton")
edata = transData(paste0(exdir,"/wee_SEAWIFS.csv"))
}

\bold{Step 2.} -- Create a polygon (\code{PBSmapping}: \code{\link[PBSmapping]{PolySet}})
  from one depth contour along the BC coast.
  See \code{\link[PBSplankton]{createDpoly}} for details.
\preformatted{
pdata = createDpoly(seePlot=FALSE)
}

\bold{Step 3.} -- Find satellite chorophyll events (gridded) in a polygon
  bounded seaward by a depth contour.
  See \code{\link[PBSplankton]{findDpoly}} for details.
\preformatted{
wee_SEAWIFS_filtered = findDpoly(edata, pdata, seePlot=FALSE,
  latBands =list(hecate=c(52,54), qcsound=c(51,52)),
  xlim=c(-133.5,-128), ylim=c(50.5,54.5), onam="wee_SEAWIFS_filtered")
}

\bold{Step 4.} -- Slice events, previously located in depth polygons, by year or
  (X,Y) to summarise Z by julian day.
  See \code{\link[PBSplankton]{sliceDpoly}} for details.
\preformatted{
wee_SEAWIFS_annual = sliceDpoly(edata=wee_SEAWIFS_filtered,
  slice="year", onam="wee_SEAWIFS_annual", seePlot=FALSE)
}

\bold{Step 5.} -- Smooth satellite chlorophyll data, calculate bloom periods,
  and estimate total chlorophyll in blooms by integrating under the curve.
  See \code{\link[PBSplankton]{calcChloro}} for details.
\preformatted{
calcChloro("wee_SEAWIFS_annual", xfld="julian", yfld="Chl",
  span=0.60, show=c(fit=TRUE,der1=FALSE,der2=FALSE), slcs.per.page=4, 
  outfile="wee_SEAWIFS_chloro", slice=list(year=1998:2001, 
  region=c("qcsound")), days.per=1)
}
}
\section{Incorporating PBSsatellite}{
  Ideally, \pkg{PBSplankton} should work with \pkg{PBSsatellite}
  (see \href{http://roam.macewan.ca/islandora/object/gm:1451}{User's Guide}).
  Preliminary steps to convert the events data in Step 1 above can be 
  achieved using \code{PBSsatellite (>=0.2.0)}, which has additional
  functionality from the original that allows the creation of \code{ncdfData}
  files from ASCII files. The process -- save R data object as a comma-
  sepearted file, convert the \code{.csv} file to a Network Common Data Format
  (NetCDF), read in the NetCDF file (\code{.nc}) to an \code{ncdfData} object
  -- is a tad inelegant, but it works using currently available tools.

  Some rudimentary steps are:
\preformatted{
require(PBSsatellite)

## Create an event data object as in Step 1 above.
## -----------------------------------------------
exdir = system.file("examples",package="PBSplankton")
edata = transData(paste0(exdir,"/wee_SEAWIFS.csv"))
write.csv(edata,"edata.csv", row.names=FALSE)

## Convert the ASCII data file to a NetCDF file and 
## import to an \code{ncdfData} object.
## -----------------------------------------------
convert.ncdfData("edata.csv", zfld="Chl", nc.filename="edata.nc", 
   summary.func=sum, offset=c(0.05,0), mv=-99, 
   dataType="Chlorophyll", dataUnits="mg/m3")
ncdata = read.ncdfData(filename="edata.nc")
print(ncdata)

## Plot a random slice of pizza
## -----------------------------------------------
pizza = sapply(ncdata,function(x){!all(is.na(x$data))})
slice = (1:length(pizza))[pizza]
plot(ncdata,slice=sample(slice,1),xlim=c(-133.5,-127.5),ylim=c(50,55))
data(nepacLL)
addPolys(nepacLL, col="gainsboro")
}
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus opus}: Institute of Ocean Sciences (IOS), Sidney BC\cr
  Last modified \code{Rd: 2018-10-18}
}
\references{
Boers, N.M., Haigh, R. and Schnute, J.T. (2004)
\href{http://www.dfo-mpo.gc.ca/Library/285693.pdf}{PBS Mapping 2: developer's guide}.
\emph{Can. Tech. Rep. Fish. Aquat. Sci.} \bold{2550}: iv + 38 pp.

Lefebvre, N.R., Boers, N., Godbout, L. and Haigh, R. (2017)
\href{http://roam.macewan.ca/islandora/object/gm:1451}{PBSsatellite 1.0: user's guide}.
\emph{Tech. Rep.} \bold{MACEWANU-CMPT-TR--2017-1}, 57 pp., 
Dept. Comp. Sci., MacEwan University, Edmonton AB.

Schnute, J.T., Boers, N.M. and Haigh, R. (2004)
\href{http://www.dfo-mpo.gc.ca/Library/285683.pdf}{PBS Mapping 2: user's guide}. 
\emph{Can. Tech. Rep. Fish. Aquat. Sci.} \bold{2549}: viii + 126 pp.

Schnute, J.T., Couture-Beil, A. and Haigh, R. (2006)
\href{http://www.dfo-mpo.gc.ca/Library/326794.pdf}{PBS Modelling 1: user's guide}. 
\emph{Can. Tech. Rep. Fish. Aquat. Sci.} \bold{2674}: viii + 114 pp.
}
\keyword{package}
\seealso{
  \code{\link[PBStools:PBStools-package]{PBStools}}
  \code{\link[PBSsatellite:PBSsatellite-package]{PBSsatellite}}
}
