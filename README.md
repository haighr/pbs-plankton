## PBSplankton: Plankton from samples to satellites ##

**PBSplankton** provides an R interface for algorithms used in biological (and fisheries) oceanography. The scope of this package is by no means comprehensive. Many of the functions provide a quick way to visualize data, and in some cases perform preliminary analyses. Though oriented to users at the Pacific Biological Station (PBS), these functions may prove useful to users at other locales. The User’s Guide is organised into sections that loosely classify the functions by theme – (1) Utility, (2) Samples, and (3) Satellites. Within each section, the functions are described alphabetically. (*User Guide incomplete at present*.)

Additional functions geared toward fisheries can be found in the package **PBStools**, upon which **PBSplankton** depends. In turn, **PBStools** depends heavily on two other R packages: **PBSmapping** (Schnute et al., 2004; Boers et al., 2004) and **PBSmodelling** (Schnute et al., 2006). **PBSplankton** has a dedicated temporary working environment called `.PBSptonEnv`. Accessor functions called `ptget`, `ptcall`, and `ptprint` enable the user to get, call, and print objects from the `.PBSptonEnv` environment, while `ptput` puts (writes) objects to `.PBSptonEnv`.

This package is still in development -- use at your own risk. 
