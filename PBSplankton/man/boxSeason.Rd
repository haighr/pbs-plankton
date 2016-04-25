\name{boxSeason}
\alias{boxSeason}
\title{
Display Seasonal Patterns Using Boxplots
}
\description{
  Summarize time-dependent data using boxplots to group by a
  specified period pattern.
}
\usage{
boxSeason(fqtName="Ex01_Sample_Info", dbName="Examples.mdb",
   type="MDB", path=getwd(), fld="S", brks="M", 
   pdf=FALSE, wmf=FALSE, ioenv=.GlobalEnv)
}
\arguments{
  \item{fqtName}{name of a file, query or table in the specified database \code{dbName}; \cr
    (can also be a local data object if \code{type="FILE"}.)}
  \item{dbName}{name of a database.}
  \item{type}{type of database (one of \code{"XLS"}, \code{"MDB"}, \code{"SQL"}).}
  \item{path}{directory path to the database.}
  \item{fld}{field name in \code{fqtName} table that contains data to summarise.}
  \item{brks}{type of period breaks; choose one of one of \code{"Q"} (quarterly),
    \code{"S"} (seasonal), \code{"B"} (bimonthly), or \code{"M"} (monthly).}
  \item{pdf}{logical: if \code{TRUE}, send figure to a postscript document file \code{.pdf}.}
  \item{wmf}{logical: if \code{TRUE}, send figure to a windows metafile \code{.wmf}.}
  \item{ioenv}{input/output environment for function input data and output results.}
}
\details{
  The data file imported from a Microsoft Access table or query called 
  \code{fqtName} must contain a field called either \code{YMD} or 
  \code{Date} that depicts the date values as strings with format \code{YYYY-MM-DD}.
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Pacific Biological Station, Fisheries and Oceans Canada, Nanaimo BC.
}
\seealso{
  \code{\link[PBSplankton]{plotDiversity}},
  \code{\link[PBSplankton]{trackComp}}
}
\examples{
local(envir=.PBSptonEnv,expr={
pbsfun=function(os=.Platform$OS.type, arch=.Platform$r_arch) {
  if (os=="windows" && arch=="i386")
    boxSeason(fld="S",brks="M",path=system.file("examples",package="PBSplankton"))
  else 
    showMessage("Only functional for Microsoft Access on 32-bit Windows OS")
  invisible() }
pbsfun()
})
}
\keyword{hplot}
\keyword{category}

