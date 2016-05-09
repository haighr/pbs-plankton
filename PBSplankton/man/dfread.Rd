\name{dfread}
\alias{dfread}
\title{
  Data Frame Read
}
\description{
  Read a data.frame using the \code{fread} function from the package \pkg{data.table}.
}
\usage{
dfread(input, ...)
}
\arguments{
  \item{input}{Either the file name to read (containing no \\n character),
    a shell command that preprocesses the file (e.g. \code{fread("grep blah filename")}) 
    or the input itself as a string (containing at least one \\n), see examples. 
    In both cases, a length 1 character string. A filename input is passed through
    \code{path.expand} for convenience and may be a URL starting http:// or file://.}
  \item{\dots}{additional arguments to pass to \code{fread}.}
}
\details{
  This function is simply a wrapper for \code{fread}.
  The only thing it does is strip the \code{"data.table"}
  class from the object.
}
\value{
  An object read by \code{fread} minus the \code{data.table} class.
}
\author{
  Rowan Haigh, Research Biologist, \cr
  Fisheries and Oceans Canada, Pacific Biological Station, Nanaimo BC
}
\seealso{
  \code{\link[data.table]{fread}},
  \code{\link[PBSplankton]{transData}}
}
\examples{
local(envir=.PBSptonEnv,expr={
  exdir = system.file("examples",package="PBSplankton")
  wiffy = dfread(paste0(exdir,"/wee_SEAWIFS.csv"))
  print(wiffy[,1:7])
})
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}

