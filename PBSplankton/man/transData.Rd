\name{transData}
\alias{transData}
\title{
  Transpose Tabulated Data from File
}
\description{
  Transpose tabulated geo-referenced data from a file into a data frame.
  The routine attempts to automatically detect the longitude and latitude fields.
}
\usage{
transData(x, zprefix="V", onames=c(x="X", y="Y", z="Chl"), byrow=TRUE)
}
\arguments{
  \item{x}{string name of the flat-table data file 
    (e.g., comma-delimited text file).}
  \item{zprefix}{string prefix before field names that essentially
    describe the date and sampling interval.}
  \item{onames}{named character vector containg the output field names
    for longitude, latitude, and z-value (e.g., chlorophyll).}
  \item{byrow}{logical, if \code{TRUE}, re-order the data.frame by the
    fields \code{cell} and \code{date}.}
}
\details{
  This function creates a long-form table from a flat table where 
  every column in the file describes gridded satellite observations
  for one time slice. The long-form table is R-friendly where every
  field describes a variable (that is, date only occurs in one column
  as opposed to hundreds of columns).
}
\value{
  Data frame with the following fields:
  \tabular{ll}{
    \code{X}        \tab Longitude (the actual field name is taken from \code{onames["x"]}.\cr
    \code{Y}        \tab Latitude (the actual field name is taken from \code{onames["y"]}.\cr
    \code{date}     \tab Date of the time slice.\cr
    \code{interval} \tab Number of days summarised in the time slice.\cr
    \code{Chl}      \tab Satellite observation of interest (the field name is taken from \code{onames["z"]}.\cr
    \code{cell}     \tab The satellite grid cell number in which an observation occurs.\cr
  }
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC
}
\seealso{
  \code{\link[PBSplankton]{dfread}},
  \code{\link[PBSplankton]{diffGeo}}
}
\examples{
local(envir=.PBSptonEnv,expr={
  exdir = system.file("examples",package="PBSplankton")
  spiffy = transData(paste0(exdir,"/wee_SEAWIFS.csv"))
  print(spiffy[1:20,])
})
}
\keyword{IO}

