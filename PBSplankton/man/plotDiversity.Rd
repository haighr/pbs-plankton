\name{plotDiversity}
\alias{plotDiversity}
\title{
  Plot Diversity of Phytoplankton Samples
}
\description{
  Plot sample diversity (\emph{e.g.}, Shanon-Wiener diversity index \eqn{H}) as a barplot 
  with an overlay of points using some other index (\emph{e.g.}, species richness \eqn{S}).
}
\usage{
plotDiversity(fqtName="Ex01_Sample_Info",dbName="Examples.mdb",
   type="MDB", path=getwd(), bars="H", pnts="S", xnames=c("SID","Batch"),
   xint=40, clrs=c("skyblue","gold","blue","green4"), addlowess=TRUE,
   f=0.2, pdf=FALSE, wmf=FALSE, ioenv=.GlobalEnv)
}
\arguments{
  \item{fqtName}{string specifying name of file, query, or table; \cr
    (can also be a local data object if \code{type="FILE"}.)}
  \item{dbName}{string specifying the name of a database.}
  \item{type}{type of database (e.g., \code{MDB} = Microsoft ACCESS database query or table).}
  \item{path}{directory path to the database.}
  \item{bars}{string value denoting field name with an index for barplots.}
  \item{pnts}{string value denoting field name with an index for points.}
  \item{xnames}{vector of strings denoting fields to use as barplot labels.}
  \item{xint}{integer value specifying number of labels to display (useful when number of bars is high).}
  \item{clrs}{vector of colours to use for (1) bars, (2) points, (3) lowess line through bars,
    (4) lowess line through points.}
  \item{addlowess}{logical: if \code{TRUE}, add lowess-fitted lines through bars and points.}
  \item{f}{numeric value for smoother span in \code{lowess}. This gives the proportion 
    of points in the plot which influence the fit at each value. Larger values give more smoothness.}
  \item{pdf}{logical: if \code{TRUE}, send figure to a postscript document file (\code{.pdf}).}
  \item{wmf}{logical: if \code{TRUE}, send figure to a windows metafile (\code{.wmf}).}
  \item{ioenv}{input/output environment for function input data and output results.}
}
\details{
  Creates a barplot of one index with overlaid points of another. The arguments 
  \code{bars} and \code{pnts} can each be one of \eqn{H} = diversity index, \eqn{S} =
  species richness, or \eqn{E} = species evenness.
}
\note{
  The Shannon-Wiener diversity index \eqn{H} is calculated
  
  \deqn{ H = \sum_{i=1}^S{p_i \log_2 p_i} }{%
         H = \Sigma p[i] * log2(p[i]);  i = 1,...,S }
  
  where \eqn{S} = species richness (number of species observed) and \cr
  \eqn{p_i}{p[i]} = proportion of each species \eqn{i} in the sample.
  
  Species evenness \eqn{E} simply measures how evenly the species are distributed
  and is calculated
  
  \deqn{ E = \frac{H}{\log_2S} }{%
         E = H / log2(S) } 

  The value of \eqn{E} ranges from 0 for situations with only a single species present 
  to 1 where all species are evenly distributed.
}
\value{
  No value is explicitly return; however, objects of interest are stored in the 
  global list object \code{PBStool}.
}
\references{
  Haigh, R. and Schnute, J.T. (2003) The Longspine Thornyhead fishery along the west coast of 
  Vancouver Island, British Columbia, Canada: portrait of a developing fishery. 
  \emph{North American Journal of Fisheries Management} \bold{23}: 120--140.
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC.
}
\seealso{
  \code{\link[PBSplankton]{boxSeason}},
  \code{\link[PBSplankton]{trackComp}}
}
\examples{
local(envir=.PBSptonEnv,expr={
pbsfun=function(os=.Platform$OS.type,arch=.Platform$r_arch) {
  if (os=="windows" && arch=="i386")
    plotDiversity(bars="H",pnts="S",xnames="Batch", 
      path=system.file("examples",package="PBSplankton"))
  else 
    showMessage("Only functional for Microsoft Access on 32-bit Windows OS")
  invisible() }
pbsfun()
})
}
\keyword{hplot}

