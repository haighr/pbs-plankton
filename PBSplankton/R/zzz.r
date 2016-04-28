# Taking cue from Roger Bivand's maptools:
.PBSptonEnv <- new.env(FALSE, parent=globalenv())  # be sure to exportPattern("^\\.PBS") in NAMESPACE

.onAttach <- function(lib, pkg)
{
	pkg_info = utils::sessionInfo( package="PBSplankton" )$otherPkgs$PBSplankton
	if( is.character( pkg_info$Packaged ) )
		pkg_date <- strsplit( pkg_info$Packaged, " " )[[1]][1]
	else
		pkg_date  <- date()

	userguide_path <- system.file( "doc/PBSplankton-UG.pdf", package = "PBSplankton" )
	year <- substring(date(),nchar(date())-3,nchar(date()))

	packageStartupMessage("
-----------------------------------------------------------
PBS Plankton ", pkg_info$Version, " -- Copyright (C) 2013-",year," Fisheries and Oceans Canada

A complete user guide 'PBSplankton-UG.pdf' is located at 
", userguide_path, "

Packaged on ", pkg_date, "
Pacific Biological Station, Nanaimo

All available PBS packages can be found at
https://github.com/pbs-software

Plankton -- from samples to satellites.
-----------------------------------------------------------

")
}
.onUnload <- function(libpath) {
	rm(.PBSptonEnv)
}

# No Visible Bindings
# ===================
if(getRversion() >= "2.15.1") utils::globalVariables(names=c(
		".findSquare",
		"addLabel", "addLegend",
		"clipVector",
		"deriFunc",
		"expandGraph",
		"hessFunc",
		"lisp",
		"packList", "pad0", "PBSdat", "predFunc",
		"resetGraph",
		"tcall", "tget", "tprint", "tput",
		"unpackList",
		"ylim", "ypos", "ytck"
	), package="PBSplankton")

