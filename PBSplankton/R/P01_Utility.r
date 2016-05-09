#===============================================================================
# Module 1: Utility
# -----------------
#  ptget...........Provide wrappers for PBSmodelling functions tget/tcall/tprint/tput/lisp
#  dfread..........Read a data.frame using data.table::fread
#  transData.......Transpose tabulated geo-referenced data from a file into a data frame.
#
#-----Supplementary hidden functions-----
#===============================================================================


#dfread---------------------------------2016-05-06
# Read a data.frame using data.table::fread
#-----------------------------------------------RH
dfread = function(input, ...){
	pout = fread(input, ...)
	attr(pout,"class") = setdiff(class(pout),"data.table")
	return(pout)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~dfread


#ttget----------------------------------2016-04-25
# Provide PBStools wrappers for PBSmodelling functions tget/tcall/tprint/tput/lisp
#-----------------------------------------------RH 
ptget   = function(...) {tget  (..., penv=parent.frame(), tenv=.PBSptonEnv)}
ptcall  = function(...) {tcall (..., penv=parent.frame(), tenv=.PBSptonEnv)}
ptprint = function(...) {tprint(..., penv=parent.frame(), tenv=.PBSptonEnv)}
ptput   = function(...) {tput  (..., penv=parent.frame(), tenv=.PBSptonEnv)}
plisp   = function(...) {lisp  (..., pos =.PBSptonEnv)}


#transData------------------------------2016-05-06
# Transpose tabulated geo-referenced data
# into a data frame; auto-detect Lon/Lat.
#-----------------------------------------------RH
transData = function(x, zprefix="V", onames=c(x="X",y="Y",z="Chl"), byrow=TRUE)
{
	on.exit(gc(verbose=FALSE))
	#x = as.character(substitute(x))
	#xdata = read.table(x,header=header,sep=sep)
	xdata = dfread(x)
	xflds = names(xdata)
	Lflds = xflds[grep("^[Ll]on|^[Ll]at|^X$|^Y$",xflds)] # location fields
	if (length(Lflds)<=1) stop("Cannot identify any Lon/Lat fields")
	if (length(Lflds)>2){
		if (length(grep("^X$|^Y$",Lflds))==2)
			lflds = Lflds[grep("^X$|^Y$",Lflds)]          ## (X,Y) fields
		else if (length(grep("^lon|^lat",Lflds))==2)
			lflds = Lflds[grep("^lon|^lat",Lflds)]        ## (lon,lat) fields
		else
			lflds = Lflds[grep("^Lon|^Lat",Lflds)]        ## (Lon,Lat) fields
	} else
		lflds = Lflds
	xfld  = lflds[grep("^[Ll]on|^X$",lflds)]
	yfld  = lflds[grep("^[Ll]at|^Y$",lflds)]
	zgrep = eval(parse(text=paste0("grep(\"^",zprefix,"\",xflds)")))
	zflds = xflds[zgrep]                                ## z-data fields (e.g., chla)
	if (length(zflds)==0) stop(paste0("Cannot identify any data fields beginning with `",zprefix,"'"))
	zdata = list(x=numeric(),y=numeric(),date=character(),interval=integer(),z=numeric(),cell=integer())
	names(zdata)[is.element(names(zdata),c("x","y","z"))] = onames
#browser();return()
	for(z in zflds) {
		zbits = strsplit(sub(zprefix,"",z),split="_")
		zdate = zbits[[1]][1]
		zdate = paste0(substring(zdate,1,4),"-",substring(zdate,5,6),"-",substring(zdate,7,8))
		zint  = as.integer(zbits[[1]][2])
		zdata[[onames["x"]]] = c(zdata[[onames["x"]]], xdata[,xfld])
		zdata[[onames["y"]]]  = c(zdata[[onames["y"]]],  xdata[,yfld])
		zdata[["date"]] = c(zdata[["date"]], rep(zdate,nrow(xdata)))
		zdata[["interval"]] = c(zdata[["interval"]],rep(zint,nrow(xdata)))
		zdata[[onames["z"]]]  = c(zdata[[onames["z"]]], xdata[,z])
		zdata[["cell"]] = c(zdata[["cell"]], 1:nrow(xdata))
	}
#browser();return()
	zdata = data.frame(zdata)
	if (byrow) zdata = zdata[order(zdata$cell,zdata$date),]
	return(zdata)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~transData


