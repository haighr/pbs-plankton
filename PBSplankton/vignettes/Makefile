## twice seems enough for now -- could also use texi2dvi, not portably
../inst/doc/PBSplankton-UG.pdf: PBSplankton-UG.Rnw
	pdflatex PBSplankton-UG.Rnw
	pdflatex PBSplankton-UG.Rnw
	mv PBSplankton-UG.pdf ../inst/doc
	@rm -f PBSplankton-UG.aux PBSplankton-UG.log PBSplankton-UG.out
