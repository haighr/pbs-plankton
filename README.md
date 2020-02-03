## PBSplankton: Plankton from samples to satellites ##
&copy; Fisheries and Oceans Canada (2013-2020)

**PBSplankton** provides an R interface for algorithms used in biological (and fisheries) oceanography. The scope of this package is by no means comprehensive. Many of the functions provide a quick way to visualize data, and in some cases perform preliminary analyses. Though oriented to users at the Pacific Biological Station (PBS), these functions may prove useful to users at other locales. The User’s Guide is organised into sections that loosely classify the functions by theme – (1) Utility, (2) Samples, and (3) Satellites. Within each section, the functions are described alphabetically. (*User Guide incomplete at present*.)

Additional functions geared toward fisheries can be found in the package **PBStools**, upon which **PBSplankton** depends. In turn, **PBStools** depends heavily on two other R packages: **PBSmapping** (Schnute et al., 2004; Boers et al., 2004) and **PBSmodelling** (Schnute et al., 2006). **PBSplankton** has a dedicated temporary working environment called `.PBSptonEnv`. Accessor functions called `ptget`, `ptcall`, and `ptprint` enable the user to get, call, and print objects from the `.PBSptonEnv` environment, while `ptput` puts (writes) objects to `.PBSptonEnv`.

**PBSplankton** represents just one of a series of R packages developed at the Pacific Biological Station (<a href="http://www.pac.dfo-mpo.gc.ca/science/facilities-installations/index-eng.html#pbs">PBS</a>) in Nanaimo, British Columbia. A more advanced version of **PBSplankton** might be available at <a href="https://github.com/pbs-software">pbs-software on GitHub</a>. Any evolving package (Windows binary and source tarball) is built after using CRAN's rigorous `R CMD check --as-cran` routine (on a Windows system) and posted to <a href="https://drive.google.com/drive/folders/0B2Bkic2Qu5LGOGx1WkRySVYxNFU?usp=sharing">Google Drive</a>. Most of the time, the revision on <a href="https://github.com/pbs-software/pbs-plankton">GitHub</a> can be built (supposedly) in R using `devtools::install_github("pbs-software/pbs-plankton")`; however, not every revision has been checked for CRAN worthiness.

This package is still in development -- use at your own risk.

As with any freely available product, there is no warranty or promise that **PBSplankton** will perform adequately for all circumstances. Additionally, coding errors are possible, and users should contact the package maintainer if bugs are detected.

Maintainer: <a href="mailto:rowan.haigh@dfo-mpo.gc.ca">Rowan Haigh</a>

<p align="right"><img src="DFOlogo_small.jpg" alt="DFO logo" style="height:30px;"></p> 

