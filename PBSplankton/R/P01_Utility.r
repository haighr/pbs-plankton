#===============================================================================
# Module 1: Utility
# -----------------
#  ptget...........Provide wrappers for PBSmodelling functions tget/tcall/tprint/tput/lisp
#
#-----Supplementary hidden functions-----
#===============================================================================

#ttget----------------------------------2016-04-25
# Provide PBStools wrappers for PBSmodelling functions tget/tcall/tprint/tput/lisp
#-----------------------------------------------RH 
ptget   = function(...) {tget  (..., penv=parent.frame(), tenv=.PBSptonEnv)}
ptcall  = function(...) {tcall (..., penv=parent.frame(), tenv=.PBSptonEnv)}
ptprint = function(...) {tprint(..., penv=parent.frame(), tenv=.PBSptonEnv)}
ptput   = function(...) {tput  (..., penv=parent.frame(), tenv=.PBSptonEnv)}
plisp   = function(...) {lisp  (..., pos =.PBSptonEnv)}


