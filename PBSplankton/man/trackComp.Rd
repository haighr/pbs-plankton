\name{trackComp}
\alias{trackComp}
\title{
  Track Composition of Phytoplankton over Time
}
\description{
  Track the mean composition of phytoplankton groups over time using 
  a specified binning interval.
}
\usage{
trackComp(fqtName=c("Ex01_Sample_Info","Ex02_Species_Abundance"),
   dbName="Examples.mdb", type="MDB", path=getwd(), ndays=15, groups=NULL, 
   dlim=c("2007-01-01",substring(Sys.time(),1,19)),
   clrs=c("green","dodgerblue","orange","yellow","red","grey"),
   pdf=FALSE, wmf=FALSE, ioenv=.GlobalEnv, ...)
}
\arguments{
  \item{fqtName}{string vector specifying table names containing sample information
    and species abundance, respectively; \cr
    (can also be two local data objects if  \code{type="FILE"}.)}
  \item{dbName}{string specifying the name of a database.}
  \item{type}{type of database (e.g., \code{MDB} = Microsoft ACCESS database query or table).}
  \item{path}{Directory path to the database.}
  \item{ndays}{integer value specifying a time interval in days to partition a year.}
  \item{groups}{vector of groups for summary; the default is: \cr
    \code{c("Skeletonema", "Thalassiosira", "Chaetoceros",
    "Other Diatoms", "Phytoflagellates", "Grazers")} }
  \item{dlim}{date limits to restrict time line of compositional changes.}
  \item{clrs}{vector of colours to shade \code{groups}.}
  \item{pdf}{logical: if \code{TRUE}, send figure to a postscript document file \code{PDF}.}
  \item{wmf}{logical: if \code{TRUE}, send figure to a windows metafile \code{WMF}.}
  \item{ioenv}{input/output environment for function input data and output results.}
  \item{\dots}{additional vectors to use for data selection 
    (e.g., \code{co="CPA", use=1}).}
}
\details{
  Creates a stacked line plot of relative composition frequencies \code{0:100}
  for each \code{group}. The polygons are coloured by \code{clrs}.
}
\value{
  No value is explicitly returned; however, objects of interest are store in the 
  global list object \code{PBStool}.
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC.
}
\seealso{
  \code{\link[PBSplankton]{boxSeason}},
  \code{\link[PBSplankton]{plotDiversity}},
  \code{\link[PBStools]{convYP}}
}
\examples{
local(envir=.PBSptonEnv,expr={
pbsfun=function(os=.Platform$OS.type,arch=.Platform$r_arch) {
  if (os=="windows" && arch=="i386")
    trackComp(ndays=15,path=system.file("examples",package="PBSplankton"))
  else 
    showMessage("Only functional for Microsoft Access on 32-bit Windows OS")
  invisible() }
pbsfun()
})
}
\keyword{hplot}

