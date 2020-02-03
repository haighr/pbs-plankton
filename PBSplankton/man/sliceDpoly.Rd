\name{sliceDpoly}
\alias{sliceDpoly}
\title{
  Slice Events in Depth Polygons
}
\description{
  Slice events, previously located in depth polygons, by 
  year or (X,Y) to summarise Z by julian day.
}
\usage{
sliceDpoly(edata, zfld="Chl", slice="year", novalue=-99,
   tooBig=30, seePlot=TRUE, onam="annual", delim="_", hadley=FALSE)
}
\arguments{
  \item{edata}{\code{list} - list of regional \code{EventData} objects (located in polygons
    by \code{findDpoly}) -- usually gridded satellite cholorophyll readings.}
  \item{zfld}{\code{character} -- field name in event data that contains chlorophyll readings.}
  \item{slice}{\code{character} -- field(s) in \code{edata} used to aggregate \code{zfld} 
    (e.g., chlorophyll) when calculating a mean for each julian day.}
  \item{novalue}{\code{numeric} -- value in \code{zfld} that indicates a missing value.}
  \item{tooBig}{\code{numeric} -- value in \code{zfld} that is unrealistically large; set to \code{NA}.}
  \item{seePlot}{\code{logical} -- if \code{TRUE}, plot the results on a map.}
  \item{onam}{\code{character} -- output name (without extension) for the binary results file (\code{*.rda}).}
  \item{delim}{\code{character} -- delimiter for file names and for concatenating field names when generating slices.}
  \item{hadley}{\code{logical} -- if \code{TRUE}, use \code{reshape}'s functions in \code{crossTab}.}
}
\details{
  Slices are often described as geo-referenced Z values (e.g., 
  chlorophyll, sea surface temperature) for a particular date (or 
  \emph{n}-day period). Here, this function redefines the slice to 
  group the Z observations over a spatio-temporal aggregate (e.g., year 
  or [X,Y]) and calculate the mean Z for each julian day represented
  in the dataset.
}
\value{
  Invisibly returns a data frame suitable as input for the function
  \code{calcChloro}. This data frame contains the fields:
  \tabular{ll}{
  \code{year}   \tab field name comes from the \code{slice} argument in the function call;\cr
  \code{X}      \tab mean X-coordinate (longitude);\cr
  \code{Y}      \tab mean Y-coordinate (latitude);\cr
  \code{julian} \tab Julian day corresponding to the original time slice in \code{edata};\cr
  \code{Chl}    \tab field name comes from the \code{zfld} argument in the function call;\cr
  \code{region} \tab automatically determines the areas from the named list object \code{edata}.
  }
  
  The function also writes binary objects \code{*.rda} to the user's
  working directory -- (i) the data frame above named using the
  \code{onam} argument, and (ii) a summary table of mean Z by Julian
  day and the user-specified \code{slice} and named \code{<onam>.summary.rda}.
}
\author{
  \href{mailto:rowan.haigh@dfo-mpo.gc.ca}{Rowan Haigh}, Program Head -- Offshore Rockfish\cr
  Pacific Biological Station (PBS), Fisheries & Oceans Canada (DFO), Nanaimo BC\cr
  \emph{locus: Institute of Ocean Sciences (IOS), Sidney BC}\cr
  Last modified \code{Rd: 2018-06-07}\cr
}
\note{
  If the slice contains some spatial element (\code{X} or \code{Y})
  then the mean Z in the results object (\code{<onam>.rda}) will 
  conatin (X,Y) coordinates weighted by Z. These coordinates are 
  currently not used for anything in particular.
}
\seealso{
  \code{\link[PBSplankton]{transData}},
  \code{\link[PBSplankton]{createDpoly}},
  \code{\link[PBSplankton]{findDpoly}},
  \code{\link[PBSplankton]{calcChloro}}
}
\examples{
local(envir=.PBSptonEnv,expr={
  ## If you have the filtered data binary file, just load it
  exdir = system.file("examples",package="PBSplankton")
  edata = transData(paste0(exdir,"/wee_SEAWIFS.csv"))
  pdata = createDpoly(seePlot=FALSE)
  wee_SEAWIFS_filtered = findDpoly(edata, pdata, seePlot=FALSE,
    latBands =list(hecate=c(52,54), qcsound=c(51,52)),
    xlim=c(-133.5,-128), ylim=c(50.5,54.5), onam="wee_SEAWIFS_filtered")
  wee_SEAWIFS_annual = sliceDpoly(edata=wee_SEAWIFS_filtered,
    slice="year", onam="wee_SEAWIFS_annual")
  view(wee_SEAWIFS_annual,10,random=TRUE)
})
}
\keyword{manip}
\keyword{hplot}

