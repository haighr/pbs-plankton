#===============================================================================
# Module 3: Satellites
# --------------------
#  calcChloro......Calculate chlorophyll peaks and integrate.
#  calcDeriv.......Calculate derivatives of loess-smoothe curve.
#  createDpoly.....Create a depth polygon from a single isobath.
#  diffGeo.........Difference a geospatial object (in km).
#  findDpoly.......Find satellite chlorophyll events in a depth polygon.
#  sliceDpoly......Slice events by year or (X,Y) to summarise Z by julian day.
#===============================================================================


#calcChloro-----------------------------2016-05-11
# Calculate chlorophyll peaks and integrate
#-----------------------------------------------RH
calcChloro  = function(dat, wd=getwd(), xfld="period", yfld="Chlo", 
  xint=0.1, span=0.75, diff.required=FALSE, slice=NULL, days.per=8,
  salmEnd="Jun-15", show=c(fit=TRUE,der1=FALSE,der2=FALSE),
  slcs.per.page=1, plts.hide=FALSE, pdf =FALSE, outfile="ChloroFits")
{
	oldpar = if (dev.cur()>1) par(no.readonly = TRUE) else NULL
	oldfac = options()$stringsAsFactors
	options(stringsAsFactors = FALSE)
	leave.idaho = function(){
		if (!is.null(oldpar)) par(oldpar)
		options(stringsAsFactors = oldfac)
		gc(verbose=FALSE)
		pdf.on = which(names(dev.list())=="pdf")
		if (any(pdf.on))
			for (i in rev(pdf.on)) dev.off(dev.list()[i])
	}
	on.exit (leave.idaho())

	## Get data
	infile = as.character(substitute(dat))
#browser();return()
	satellite = strsplit(infile,"\\.|_|-|\\||\\+")[[1]][1]  ## assumes leading string before first "." is satellite name
	if (!is.data.frame(dat)) {
		if (is.character(dat) && length(dat)==1){
			expr = paste0("getFile(\"",dat,"\",path=wd,reload=TRUE,tenv=penv()); dat=",dat)
			eval(parse(text=expr))
		} else
			stop ("`dat' must be a data frame or the name of a data file")
	}
	rawData = dat # save just in case
	if (any(grepl("[lL]on",names(dat)))) {
		dat$sign.lon = sign(dat[,grep("[lL]on",names(dat))[1]])
		dat[,grep("[lL]on",names(dat))[1]] = abs(dat[,grep("[lL]on",names(dat))[1]])
	}
	## Get slice
	if (is.null(slice))
		dat$slice = "use-all-records"
	else if ("slice" %in% names(dat) && is.vector(slice) && is.character(slice)) {
		dat = dat[is.element(dat$slice,slice),]
		sflds = "slice"
	} else {
		##Extract the slice
		if (is.list(slice)) {
			##Grab the specified slice
			for (i in 1:length(slice)) {
				ii=names(slice)[i]; iii=slice[[i]]
				dat = dat[is.element(dat[[ii]],iii),]
				if (nrow(dat)==0) {showMessage("No data for specified slice"); return("No data") }
			}
			sflds = names(slice)
		}
		else {
			##Create a slice if field names are given instead of a list
			flds = names(dat)
			sflds = slice[slice %in% flds]
		}
		dat$slice = apply(dat[,sflds],1,paste,collapse="-")
	}
	slice = sort(unique(dat$slice))
	
	res.slice = list()
	stanzas = list()

	rc = .findSquare(slcs.per.page); rc[1] = rc[1]*sum(show)
	if (pdf) {
		## http://stackoverflow.com/questions/14366406/determine-ghostscript-version
		gsversion <- names(which(Sys.which(c("gs", "gswin32c", "gswin64c")) != ""))
		if (length(gsversion)==0) stop("Install Ghostscript and make sure the bin directory is on the path")
		else gsverion = gsversion[1]

		pdf(file=paste0(outfile,".pdf"),onefile=TRUE,width=8.5,height=11,paper="special")
		#pdf(file=paste0(outfile,".pdf"),onefile=TRUE,width=6,height=switch(sum(show),5,8,9),paper="special")
		## gswin64c has a command line limit so create bookmarks in stanzas of 50; also need to call backwards on command line
		## start building a batch file to bookmark the PDF (https://groups.google.com/forum/#!topic/comp.text.pdf/TslRCZH6x70)
		bookhead = c(
			#"%R_GSCMD%\\gswin64c.exe ^",
			paste0(gsversion, " ^"),
			paste0(" -o outfile.pdf ^"),
			" -sDEVICE=pdfwrite ^")
		bookend =  " -f infile.pdf"
	}
	if (sum(show)==1) par(mfrow=rc)
	else              par(mfcol=rc)
	par(mar=c(3,3,2,0.5), oma=c(0.5,0.5,0.5,0.5), mgp=c(1.6,0.5,0), cex=ifelse(pdf,0.8,1))
	sss = as.character() # initialize title text for one summary title (not used currently)
#browser();return()

	## Start the output file (comma delimited)
	catfile = paste0(outfile,".csv")
	cat(paste0("File = '",infile,"' -- Summary of fits using loess with span=",span,"\n"),file=catfile)
	#collect = c(sflds,"peakX","peakY","peakX0","peakX1","peakDur","peakInt","bloomN","bloomX0","bloomX1","bloomDur","bloomInt","totalX0","totalX1","totalInt")
	collect = c(sflds,"peakX","peakY","peakX0","peakX1","peakDur","peakInt","salmInt","totalX0","totalX1","totalInt")
	cat(paste0(collect,collapse=","),"\n",sep="",file=catfile,append=TRUE)

	for (s in 1:length(slice)){ # 563 slices
		if (s%%slcs.per.page == 0) ss  = as.character()
		if (s%%slcs.per.page == 1) sss = as.character()
		ss = slice[s]
		sss = c(sss,ss)
#print(c(s,sss))#; print(c(s,ss))
		sz = ceiling(s/50)   # stanza index
		if (pdf && s%%50 == 1) bookmark = bookhead

		## Slice data
		sdat = dat[is.element(dat$slice,ss),]
		if (any(grepl("[Dd]ate",names(sdat)))) {
			sdate = sdat[,grep("[Dd]ate",names(sdat))[1]]
			days.per = as.numeric(diff(sdate)[1])
			sper0 = seq(sdate[1], as.Date(paste0(substring(sdate[1],1,4),"-01-01")), -days.per)
			sdat$period = seq(1,nrow(sdat),1) + length(sper0) -1
			xfld = "period"
		}
		sdat = sdat[!is.na(sdat[,yfld]),]
		x  = sdat[,xfld]
		y  = sdat[,yfld]
		if (diff.required)
			y = c(y[1],diff(y))
		nx = length(x)
		ymax = max(y,na.rm=TRUE)
		cutoff = 1.05 * median(y) # long-term median

		## Generate the statistics and the functions for the chlo-function and its derivative
		resChlo = calcDeriv(x,y,span=span)
		unpackList(resChlo)

		## Predicted values from smoothed function
		predVal = predFunc(x)
		xsmoo   = seq(x[1],rev(x)[1],xint)
		ysmoo = predFunc(xsmoo)

		## Derivative.  grad crashes near bounds
		xderi   = x[2:(nx-1)]
		xderi   = xsmoo[2:(length(xsmoo)-1)]
		deriVal = deriFunc(xderi)
		hessVal = hessFunc(xderi)
		yIntTot = stats::integrate(predFunc,min(xsmoo),max(xsmoo))$value

		## Determine peaks in smoothed trend
		dsign = sign(deriVal)
		dbrks = diff(c(dsign[1],dsign))
		#dbrks = diff(c(dsign,rev(dsign)[1]))
		xpos  = (1:length(deriVal))[dbrks==-2]
		xlook = ceiling(length(xsmoo)/(2*nx))
		zL = sapply(xpos,function(x){
			i = max(1,(x-xlook)):min((x-1),length(dbrks))
			all(dbrks[i]==0)})
		zR = sapply(xpos,function(x){
			i = (x+1):min((x+xlook),length(dbrks))
			all(dbrks[i]==0)})
		xpos = xpos[zL&zR]
		xpos = sapply(xpos,function(x){ # check y-values around inflection for the highest y
			xx = x+c(-1,0,1)
			yy = ysmoo[2:(length(ysmoo)-1)][xx]
			z  = grep(max(yy),yy)[1]
			xx[z]
			})
		xpeak = round(xderi[xpos],5) #-0.5

		## Check for continuous declines at start or continuous increase at end of series
		ysign = sign(ysmoo-cutoff)
		if (any(ysign==1)) {
			names(ysign) = xsmoo
			xsign = which(ysign > 0, useNames=TRUE)
			#http://stackoverflow.com/questions/18508363/split-a-numeric-vector-into-continuous-chunks-in-r
			xbloc = split(xsign, cumsum(seq_along(xsign) %in% (which(diff(xsign)>1)+1)))
			xbloc = sapply(xbloc,function(x){round(as.numeric(names(x)),5)},simplify=FALSE)
			if (sapply(xbloc[1],function(x,y){!any(y %in% x)},y=xpeak))
				xpeak = c(xbloc[[1]][1],xpeak)
			if (sapply(rev(xbloc)[1],function(x,y){!any(y %in% x)},y=xpeak))
				xpeak = c(xpeak,rev(rev(xbloc)[[1]])[1])
		}
		ypeak = predFunc(xpeak)

		segx = function(x){as.vector(rbind(x,x,rep(NA,length(x))))}
		segy = function(x){as.vector(rbind(rep(0.025*ymax,length(x)),x,rep(NA,length(x))))}
		olga = function(x,y,...){polygon(c(x[1],x,rev(x)[1]),(c(0,y,0)),...)}

		## Spring bloom estimation
		getBloom = function(x,y,xper,ycut){
			nper = length(x)
			i1 = 1:(nper-xper+1)
			i2 = xper:nper
			z  = apply(cbind(i1,i2),1,function(i){
				y[i[1]]:y[i[2]] > ycut})
			z = c(z,rep(rev(z)[1],xper-1))
			if (!any(z)) return(FALSE) # no bloom
			x0=c(which(z)[1],which(z)[c(1,diff(which(z)))>1])
			x9=c(which(z)[c(diff(which(z),1))>1],rev(which(z))[1])
			tlims = as.data.frame(array(x[rbind(x0,x9)],dim=c(2,length(x0)),dimnames=list(c("start","end"),names=paste0("peak",1:length(x0)))))
			blims = as.data.frame(array(y[rbind(x0,x9)],dim=c(2,length(x0)),dimnames=list(c("start","end"),names=paste0("peak",1:length(x0)))))
			peaks = sapply(1:length(x0),function(i){
				xx = x[x0[i]:x9[i]]
				yy = y[x0[i]:x9[i]]
				z  = match(max(yy),yy)[1]
				maxy = yy[z]; names(maxy) = xx[z]
				return(maxy) } )
			return(list(tlims=tlims,blims=blims,peaks=peaks))
		}
		bloom = getBloom(x=xsmoo,y=ysmoo,xper=3,ycut=cutoff)
#if(substring(ss,1,4)=="2003") {browser();return()}
		if (is.logical(bloom) && !bloom) {
			nbloo = 0
		}
		else {
			peaky    = sapply(bloom$tlims,function(i){any(xpeak>=i[1] & xpeak<=i[2])})
			nbloo    = sum(peaky)
			xbloo    = sapply(bloom$tlims[peaky],function(i){segx(i)})
			ybloo    = sapply(bloom$blims[peaky],function(i){segy(i)})
			zbloo    = list(x=bloom$tlims[peaky], y=bloom$blims[peaky])
			yIntBloo = sapply(zbloo$x,function(x){stats::integrate(predFunc,x[1],x[2])$value})
			xPkDur   = apply(zbloo$x,2,diff)
			peaks    = bloom$peaks
			mbloo    = is.element(round(peaks,5),round(max(peaks),5))
			#peaks = bloom$peaks
			#mbloo = round(ypeak[peaky],5)==round(max(ypeak[peaky]),5) # index of maximum bloom
			#maxPk = ypeak[peaky][mbloo]
		}
		## Assume salmon end date exists
		salmX1 = strptime(paste0(substring(ss,1,4),"-",salmEnd),"%Y-%b-%d")$yday+1
		if (zbloo$x[1,1] <= salmX1) { ## e.g., Jun 15 = 166 or 167
			#salmon = sapply(zbloo,function(x){matrix(x[,1],ncol=1)},simplify=FALSE)
			salmon = sapply(zbloo,function(x){data.frame(salmon=x[,1])},simplify=FALSE)
			salmon$x[2,1] = salmX1
			salmon$x[2,1] = min(salmon$x[2,1],zbloo$x[2,1])
			salmon$y[2,1] = predFunc(salmX1)
			salmBloo = sapply(salmon$x,function(x){stats::integrate(predFunc,x[1],x[2])$value})
			#salmDur = salmX1 - zbloo$x[mbloo][1,1]
			salmDur = salmon$x[2,1]-salmon$x[1,1]
			#salmon = sapply(salmon,function(x){dimnames(x)=list(c("start","end"),"salmon");x},simplify=FALSE)
			#salmBloo = stats::integrate(predFunc,salmon$x[1],salmon$x[2])$value
#browser();return()
		} else
			salmon = NULL
	
		## Plot the results
		#par(mar=c(3,3,0.5,0.5), oma=c(0.5,0.5,2,0.5), mgp=c(1.6,0.5,0), cex=ifelse(pdf,0.8,1))
		par(cex=ifelse(rc[2]>1,.8,1))
		xlim = c(1,ceiling(365/days.per)); ylim=c(0,max(1.10*ymax,1.75*cutoff))

		### Continue building the bookmarks
		titext = if (!exists("sflds")) ss else paste0(sflds,"=", strsplit(ss,split="-")[[1]], ", ",collapse="")
		titext = substring(titext,1,nchar(titext)-2) # get rid of trailing ", "
		titext = paste0(toupper(satellite),": ",titext)
		if (pdf && s%%slcs.per.page == 0) {
			#sssmat = sapply(sapply(sss,strsplit,split="-"),function(x){return(as.numeric(x))})
#if (s==12) {browser();return()}
			sssmat = sapply(sapply(sss,strsplit,split="-"),function(x){return(x)})
			bootext = paste0(sflds,"=",apply(sssmat,1,function(x){paste0(unique(x),collapse=",")}),"; ",collapse="")
			bootext = substring(bootext,1,nchar(bootext)-2) # get rid of trailing "; "
		#if (pdf){
			bookmark = c(bookmark,
			paste0(" -c \"[/Page ",round(s/slcs.per.page)," /View [/XYZ null null null] /Title (",bootext,") /OUT pdfmark\" ^"))
			if (s%%(50*slcs.per.page) == 0){
				bookmark = c(bookmark, bookend)
				stanzas[[sz]] = bookmark
				#bookmark = c(bookmark, gsub("outfile",paste0("temp",stanza),bookend))
				#bookmark = c(bookmark, gsub("01\\\\.pdf",paste0(stanzas[sz+1],"\\.pdf",bookhead)))
			}
		}

		if (show["fit"]) {
			## Plot original data and smoothed function
			plot(x,y,xlim=xlim,ylim=ylim,xlab=xfld,ylab=yfld,type="n")
			#abline(v=167,lty=2,col="lightgrey")
			#points(salmX1,par()$usr[3]-diff(par()$usr[3:4])*0.02,pch=24,bg="salmon",cex=1,xpd=NA)

			if (nbloo>0) {
				for (i in 1:nbloo){ 
					xpoly = seq(zbloo$x[1,i],zbloo$x[2,i],xint)
					ypoly = ysmoo[is.element(round(xsmoo,5),round(xpoly,5))]
					olga(xpoly,ypoly,col="honeydew",border=FALSE)
				}
				for (i in 1:nbloo){ 
					lines(xbloo[,i],ybloo[,i],col="deepskyblue",lty=2,lwd=2)
					text(zbloo$x[,i],rep(0,2),labels=round(zbloo$x[,i],1),col="blue",adj=c(0.5,1),cex=0.8)
				}
			}
			abline(h=cutoff,col="gainsboro",lwd=2)
			points(x,y,pch=20,col="red",cex=1.5)
			lines(segx(xpeak),segy(ypeak),col="dodgerblue",lwd=2)
			text(xpeak,rep(0,length(xpeak)),labels=round(xpeak,1),col="blue",adj=c(0.5,1),cex=0.8)
			lines(xsmoo,ysmoo,col="green4",lwd=2)
			legtxt = paste0("Fitted trend (span=",span,"):  ")
			legtxt = paste0(legtxt, ifelse(rc[2]==1,"\n",""),
				paste0("\u03A3year = (",round(diff(range(sdat[,xfld]))+1,0),"d, ",round(yIntTot,0), "mg/m\u00B3)") )
				#paste0("\u03A3 Tchl = ",round(yIntTot,1), "\u03A3 Tdur = ", round(diff(range(sdat[,xfld]))+1,1)) )
			if (nbloo>0) {
				legtxt = paste0(legtxt, ifelse(rc[2]>1,"\n",";  "),
				"\u03A3bloom = (", round(sum(xPkDur)*days.per,0),"d, ",round(sum(yIntBloo),0), "mg/m\u00B3)")
				#"\u03A3 Bchl = ", round(sum(yIntBloo),1),
				#"; \u03A3 Bdur = ", round(sum(xPkDur),1), ifelse(days.per==1, " days",paste0(" periods (", round(sum(xPkDur)*days.per,1), " days)")) )
			}
			if (!is.null(salmon)) {
				legtxt = paste0(legtxt, ifelse(rc[2]>=1,";  ","\n"),
				"\u03A3salmon = (", round(sum(salmDur)*days.per,0),"d, ",round(sum(salmBloo),0), "mg/m\u00B3)")
				#ifelse(rc[2]>1,"\n","; "),"\u03A3 Schl = ", round(sum(salmBloo),1),
				#"; \u03A3 Sdur = ", round(sum(salmDur),1), ifelse(days.per==1, " days",paste0(" periods (", round(sum(salmDur)*days.per,1), " days)")) )
			}
			addLegend(ifelse(pdf,0.005,0.01),ifelse(rc[2]>1,0.92,0.95),lwd=2,col="green4",
				legend=legtxt,bty="n",adj=c(0,ifelse(rc[2]>1,0.75,0.75)),yjust=0,cex=ifelse(pdf,0.8,1))
			mtext(titext,side=3,line=0.25,cex=1.2,col="blue")
			box()
			polygon(x=as.vector(apply(salmon$x,1,rep,2)), y=par()$usr[3] + diff(par()$usr[3:4])*0.0075*c(-1,1,1,-1), col="salmon",border="gainsboro", xpd=NA)
		}
		if (show["der1"]) {
			##plot first derivative of smoothed function
			plot(xderi,deriVal,type="h",xlim=xlim,xlab=xfld,col="slategray1",lwd=2)
			abline(h=0,col="gainsboro",lwd=2)
			points(xpeak,rep(0,length(xpeak)),pch=21,cex=1.2,col="red",bg="gold")
			addLabel(0.95,0.9,"First derivative",cex=1.2,adj=1)
			box()
		}
		if (show["der2"]) {
			##plot second derivative of smoothed function
			plot(xderi,hessVal,type="h",xlim=xlim,xlab=xfld,col="salmon",lwd=2)
			addLabel(0.1,0.9,"Second derivative",cex=1.2,adj=0)
			abline(h=0,col="gainsboro")
			box()
		}
		outlist = list(resChlo=resChlo,predVal=predVal,deriVal=deriVal,xpeak=xpeak,ypeak=ypeak)
		if (nbloo>0) 
			outlist = c(outlist, list(bloom=bloom, peaky=peaky, zbloo=zbloo,
			yIntTot=yIntTot, yIntBloo=yIntBloo, xPkDur=xPkDur))
		res.slice[[ss]] = outlist
	
		## Dump results into outfile
		#  c(sflds,"peakY","peakX0","peakX1","peakDur","peakInt","bloomN","bloomX0","bloomX1","bloomDur","bloomInt","totalX0","totalX1","totalInt")
		sres = sapply(sflds,function(x){unique(sdat[[x]])})
		pres = if (nbloo==0)
			c(peakX="", peakY="", peakX0="", peakX1="", peakDdur="", peakInt="", salmInt="")
				#bloomN=0, bloomX0="", bloomX1="", bloomDdur="", bloomInt=0 )
		else c(
			peakX     = as.numeric(names(peaks)[mbloo]),
			peakY     = signif(as.vector(peaks[mbloo]),4), 
			peakX0    = zbloo$x[mbloo][1,], 
			peakX1    = zbloo$x[mbloo][2,],
			peakDdur  = as.vector(xPkDur[mbloo]), 
			peakInt   = as.vector(signif(yIntBloo[mbloo],4)),
			#bloomN    = nbloo, 
			#bloomX0   = min(zbloo$x["start",]), 
			#bloomX1   = max(zbloo$x["end",]), 
			#bloomDdur = sum(xPkDur),
			#bloomInt  = signif(sum(yIntBloo),4))
			salmInt    = ifelse(is.null(salmon),"",signif(salmBloo,4))
			)
		sres = c(sres, pres, 
			totalX0   = min(sdat[,xfld]),
			totalX1   = max(sdat[,xfld]),
			totalInt  = signif(yIntTot,4)
		)
		cat(paste0(sres,collapse=","),"\n",sep="",file=catfile,append=TRUE)
		if (s%%slcs.per.page == 999) {
			sssmat = sapply(sapply(sss,strsplit,split="-"),function(x){return(as.numeric(x))})
			titext = paste0(sflds,"=",apply(sssmat,1,function(x){paste0(unique(x),collapse=",")}),"; ",collapse="")
			titext = substring(titext,1,nchar(titext)-2) # get rid of trailing ", "
			mtext(titext,side=3,outer=TRUE,line=0,cex=1.2,col="blue")
		}
	}
	if (pdf){
		dev.off()
		#bookmark = c(bookmark, paste0(" -f .\\",outfile,".pdf"))
		if(s%%50 != 0) bookmark = c(bookmark,bookend)
		stanzas[[sz]] = bookmark
		names(stanzas) = paste0("stanza",pad0(length(stanzas):1,2))
		stanzas[[1]][2] =  gsub("outfile",paste0(outfile,"BM"),stanzas[[1]][2])
		lastF = grep("infile",stanzas[[length(stanzas)]])
		stanzas[[length(stanzas)]][lastF] = gsub("infile",outfile,stanzas[[length(stanzas)]][lastF])
		for (i in 1:length(stanzas)) {
			ii = names(stanzas)[i]
			iii = c(names(stanzas),"dummy")[i+1]
			if (any(grepl("outfile",stanzas[[i]]))){
				outL = grep("outfile",stanzas[[i]])
				stanzas[[i]][outL] = gsub("outfile",ii,stanzas[[i]][outL])
			}
			if (any(grepl("infile",stanzas[[i]]))){
				inL = grep("infile",stanzas[[i]])
				stanzas[[i]][inL] = gsub("infile",iii,stanzas[[i]][inL])
			}
		}
		stanzas = stanzas[length(stanzas):1]
		outlines = character()
		for (i in names(stanzas))
			outlines = c(outlines,stanzas[[i]])
		cmd = paste0(tempdir(),"/addBM.bat")
		writeLines(outlines,con=cmd)
		system2(cmd) # really slow when no.slices > 50
	}
	return(invisible(res.slice))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcChloro


#calcDeriv------------------------------2015-04-01
# Calculate derivatives (with help from Wayne Hajas)
#--------------------------------------------WH/RH
calcDeriv <- function(x,y,span=0.5)
{
	## loess to give a smooth approximation
	fitline = loess(y~x,span=span)

	## Encapsulate that within a function
	predFunc  = function(tt){
		result = predict(fitline,tt)  ## stats
		return(result)}

	## Function for the first derivative
	deriFunc  = function(tt){
		result = grad(predFunc,tt,  ## numDeriv
			method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE) )
		return(result)}

	## Function for the second derivative
	hessFunc = function(tt){
		result = grad(deriFunc,tt,  ## numDeriv
			method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE) )
		#result = numDeriv::hessian(predFunc,tt,method.args=list(d=0.1))
			#method.args=list(eps=1e-4, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE) )
		return(result)}

	## Calculate the root. grad crashes near the bounds
	#root = stats::uniroot(deriFunc,c(min(x)+.1,-.1+max(x)))#$root
	#ymax = predFunc(root$root)
	integratey = integrate(predFunc,min(x),max(x))$value  ## stats

	## Put results in a list
	result=list(
		fitline=fitline, predFunc=predFunc, deriFunc=deriFunc, 
			hessFunc=hessFunc, integratey=integratey)
	return(result)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~calcDeriv


#createDpoly----------------------------2016-05-06
# Create a depth polygon from a single isobath.
#-----------------------------------------------RH
createDpoly = function(isob=1000, corner=c("TR","BR"),
   xbox=c(-134,-134,-124.5,-124.5), ybox=c(48,54.4,54.4,48),
   seePlot=FALSE)
{
	if (!isob %in% seq(100,1800,100)) stop("Choose one iobath from 100 to 1800 by incremenst of 100")
	abox = as.PolySet(data.frame(PID=rep(1,4),POS=1:4,X=xbox,Y=ybox),projection="LL")
	amid = calcCentroid(abox)
	data(isobath, package="PBSdata", envir=penv())
	isoLine = isoPoly = isobath[is.element(isobath$PID,isob),]
	isoPoly$Y[1] = min(ybox)
	isoPoly$Y[nrow(isoPoly)] = max(ybox)
	for (i in corner) {
		z = switch(match(i,c("BL","TL","TR","BR")),
			{xbox<amid$X & ybox<amid$Y},
			{xbox<amid$X & ybox>amid$Y},
			{xbox>amid$X & ybox>amid$Y},
			{xbox>amid$X & ybox<amid$Y} )
		isoPoly = rbind(isoPoly,data.frame(c(isoPoly[nrow(isoPoly),1:3]+c(0,0,1),X=xbox[z],Y=ybox[z])) )
	}
	save("isoPoly",file="isoPoly.rda")

	if (seePlot) {
		exlim = extendrange(isoPoly$X)
		eylim = extendrange(isoPoly$Y)
		data(nepacLLhigh, package="PBSmapping", envir=penv())
		bccoast = clipPolys(nepacLLhigh,xlim=exlim,ylim=eylim)
		expandGraph()
		plotMap(bccoast,xlim=exlim,ylim=eylim,col="moccasin",plt=NULL)
		addPolys(isoPoly,border="blue",col=lucent("blue",0.2))
	}
	return(invisible(isoPoly))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~createDpoly


#diffGeo--------------------------------2016-05-05
# Transform a geographic object into a coordinate
# system with units in kilometres and an origin
# near zero or centered on the mean.
#--------------------------------------------TD/RH
diffGeo <- function(dat, offset.from.zero=0.1, 
   xfld="X", yfld="Y", center=TRUE, replace=FALSE, centroid=NULL)
{
	## constants
	dat$X = dat[,xfld]; dat$Y = dat[,yfld]
	user.centroid = centroid
	centroid =list()
	centroid[["usr"]] = user.centroid
	centroid[["old"]] = c(mean(dat$X, na.rm=TRUE), mean(dat$Y, na.rm=TRUE))
	minX  = min(dat$X) - offset.from.zero ## minimum LON
	minY  = min(dat$Y) - offset.from.zero ## minimum LAT
	
	proj  = attributes(dat)$projection
	is.deg = !is.null(proj) && proj=="LL"
	if (is.deg) {
		deg2rad =  180/pi          ## degrees to radians (if degrees)
		rminX   = minX/deg2rad
		rminY   = minY/deg2rad
		R       = 6371.2           ## radius (km) of Earth

		## geometric calculations (also see function `calcGCdist')
		convRX = function(radX,rminX,radY,R) { ## convert radians to X-values (km)
			difX   = abs(radX-rminX)
			betaX  = sin(radY)*sin(radY) + cos(radY)*cos(radY)*cos(difX)  ## hold Y constant
			betaX  = acos(betaX)
			Xe = R * betaX
		}
		convRY = function(radY,rminY,R) { ## convert radians to Y-values (km)
			betaY  = sin(radY)*sin(rminY) + cos(radY)*cos(rminY)
			betaY  = acos(betaY)
			Yn = R * betaY
		}
		radX   = dat$X/deg2rad
		radY   = dat$Y/deg2rad
		dat$Xe = convRX(radX,rminX,radY,R)
		dat$Yn = convRY(radY,rminY,R)
	} else {
		dat$Xe = dat$X
		dat$Yn = dat$Y
	}
	if (center) {
		if (!is.null(centroid[["usr"]])){
			if (is.deg) {
				uX = convRX(centroid$usr[1]/deg2rad,rminX,centroid$usr[2]/deg2rad,R)
				uY = convRY(centroid$usr[2]/deg2rad,rminY,R)
			} else {
				uX = centroid$usr[1]  ## be sure to standardise the zone across multiple
				uY = centroid$usr[2]  ## objects when converting to UTM using `convUL'.
			}
			centroid[["new"]] = c(uX, uY)
		} else
			centroid[["new"]] = c(mean(dat$Xe, na.rm=TRUE), mean(dat$Yn, na.rm=TRUE))
		dat$cXe = dat$Xe - centroid[["new"]][1]
		dat$cYn = dat$Yn - centroid[["new"]][2]
		attr(dat, "centroid") = centroid
	}
	if (replace) { ## set the transformed coordinates to the primary X & Y
		dat$X.orig = dat$X
		dat$Y.orig = dat$Y
		dat$X = if (center) dat$cXe else dat$Xe
		dat$Y = if (center) dat$cYn else dat$Yn
	}
	attr(dat,"projection") = 1 
	attr(dat, "zone") = NULL
	return(dat)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~diffGeo


#findDpoly------------------------------2016-05-11
#  Find satellite chorophyll events (gridded) in
#  a polygon bounded seaward by a depth contour.
#-----------------------------------------------RH
findDpoly = function(edata, pdata, zfld="Chl", novalue=-99, tooBig=30,
   latBands, seePlot=TRUE, xlim, ylim, onam="filtered")
{
	on.exit(gc(verbose=FALSE))

	### Modify the data
	eflds = names(edata)
	if (!all(is.element(c("X","Y"),eflds))){
		names(edata)[grep("lon",names(edata))]="X"
		names(edata)[grep("lat",names(edata))]="Y"
	}
	### Locate data in isobath polygon
	events = as.EventData(data.frame(EID=1:nrow(edata),edata),projection="LL")
	events[,zfld][is.element(events[,zfld],novalue)] = NA
	events[,zfld][events[,zfld]>tooBig & !is.na(events[,zfld])] = NA
	locset = findPolys(events, pdata, maxRows = 1e+06)
	if (length(locset)==0) stop("No events occur in the target polygon")
	EinP   = events[is.element(events$EID,sort(unique(locset$EID))),]
	if (!is.element("id",names(EinP)))
		EinP$id = paste0(EinP$Y,EinP$X)
	bucket = split(EinP[,zfld],EinP$id)
	goodchl= sapply(bucket,function(x){any(x>=0 & x<=30,na.rm=TRUE)})  ## valid chl data

	### Isolate the data into exact polygonal areas
	egood  = EinP[is.element(EinP$id,names(goodchl)[goodchl]),]
#browser();return()

	if (!missing(latBands) && !is.null(latBands)){
		bgood = list()
		for (i in 1:length(latBands)) {
			ii = names(latBands)[i]
			if (is.null(ii)) ii = paste0("band",i)
			iii = latBands[[i]]
			bhave = egood$Y<iii[2] & egood$Y>=iii[1]
			if (!any(bhave)) next
			bgood[[ii]] = egood[bhave,]
		}
	}
	if (length(bgood)==0) rm(bgood, envir=penv())

	### Plot the data positions
	if (seePlot) {
		if (missing(xlim)) xlim = extendrange(events$X)
		if (missing(ylim)) ylim = extendrange(events$Y)
		data(nepacLLhigh, package="PBSmapping", envir=penv())
		expandGraph(cex=1.5)
		plotMap(pdata,border="green",col="honeydew",xlim=xlim,ylim=ylim,plt=NULL)
		abline(h=c(51,52,54))
		addPolys(nepacLLhigh,col="moccasin")
		addPoints(EinP,pch=20,col="red")
		addPoints(egood,pch=20,col="green4")
		if (exists("bgood",envir=penv())) {
			bcol = c("blue","gold","purple","orange")
			bcol = rep(bcol,length(bgood))[1:length(bgood)]
#browser();return()
			for (i in 1:length(bgood))
				addPoints(bgood[[i]], pch=20, col=bcol[i])
		}
	}
	### Save the data
	if (exists("bgood",envir=penv()))
		ogood = bgood
	else
		ogood = list(all=egood)
	omess = paste0(onam,"=ogood; ",
		"save(\"",onam,"\",file=\"",onam,".rda\")")
	eval(parse(text=omess))
	return(invisible(ogood))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~findDpoly


#sliceDpoly-----------------------------2016-05-11
#  Slice events located in depth polygons by 
#  year or (X,Y) to summarise Z by julian day.
#-----------------------------------------------RH
sliceDpoly = function(edata, zfld="Chl", slice="year", novalue=-99,
   tooBig=30, seePlot=TRUE, onam="annual", delim="_")
{
	on.exit(gc(verbose=FALSE))

	if (!is.list(edata) || is.data.frame(edata))
		stop("User must pass a list of data.frames through the argument `edata' (see function `findDpoly')")
	anames = names(edata) ## area names

	## Summarize the data by area and year (if annual slices) or by julian day slices
	orchard=list()        ## to collect data for `calcChloro'
	crumble = list()      ## save summary tables
	zmax = 0
	for (i in anames) {
		idat = edata[[i]]
		idat[,zfld][is.element(idat[,zfld],novalue)] = NA
		idat[,zfld][idat[,zfld]>tooBig & !is.na(idat[,zfld])] = NA
		idat = idat[order(idat$date),]  ## ensure the dates are ascending
		idat$year = as.numeric(substring(idat$date,1,4))
		yrs  = .su(idat$year)
		idat$slice = apply(idat,1,function(x,y){paste0(x[y],collapse=delim)},y=is.element(names(idat),slice))
		#per  = as.numeric(names(rev(sort(table(diff(as.Date(.su(idat$date))))))[1]))       ## most common period between time slices
		#MD0  = names(table(sapply(split(idat$date,idat$year),function(x){substring(x[1],6)})))[1] ## most common starting Month-Day
		moday = .su(unlist(sapply(sapply(split(idat$date,idat$year),function(x){substring(x,6)}),.su)))
		duniq = sort(unique(moday))
		zdiff = diff(strptime(duniq,"%m-%d")$yday+1)
		zdiff = zdiff = c(zdiff[1],zdiff)
		zgood = .su(zdiff[zdiff>=(max(zdiff)-1)])
		zbad  = .su(zdiff[zdiff<(max(zdiff)-1)])
		dgood = duniq[zdiff %in% zgood]; names(dgood) = dgood
		if (length(zbad) > 0) {
		dbad  = c(duniq[1],duniq)[zdiff %in% zbad]; names(dbad)  = duniq[zdiff %in% zbad] ## index off days with on days
		} else dbad=NULL
		duse  = c(dgood, dbad)
		goday = duse[substring(idat$date,6)]
		##-----------------------------------------
		idat$julian = strptime(goday,"%m-%d")$yday+1
		uI = sort(unique(idat$slice));  nI = length(uI)
		uJ = sort(unique(idat$julian)); nJ = length(uJ)
		uX = sort(unique(idat$X));      nX = length(uX)
		uY = sort(unique(idat$Y));      nY = length(uY)
		isumm  = array(NA,dim=c(nJ,nI),dimnames=list(julian=uJ,slice=uI))
		icells = split(idat,idat$slice)
		lenv = sys.frame(sys.nframe())
		fornow = sapply(icells, function(jdat){
			zlist = split(jdat[,zfld],jdat$julian)
			zmean = sapply(zlist,function(x){if (all(is.na(x))) NA else mean(x,na.rm=TRUE)})
			ii = as.character(sort(unique(jdat$slice)))
			jj = names(zmean)
			zenv = sys.frame(sys.nframe())
			tget(isumm,penv=lenv,tenv=zenv)
			isumm[jj,ii] = zmean
			tput(isumm,penv=zenv,tenv=lenv)
		})
		zmean = apply(isumm,1,mean,na.rm=TRUE)
		zmean[is.nan(zmean)] = NA
		zmax  = max(c(zmax,zmean),na.rm=TRUE)
#browser();return()
		attr(isumm,paste0(zfld,".mean")) = zmean
		crumble[[i]] = isumm

		## Calculate (X,Y) locations weighted by julian-day cholorophyll -- only makes a difference if slice = c("X","Y")
		ctab = convCT(crossTab(idat,c("julian","slice"),"Chl",mean,na.rm=TRUE))
		XYwtd = list()
		for (xy in c("X","Y")){
			xtab = convCT(crossTab(idat,c("julian","slice"),xy,mean,na.rm=T))
			wtab = t(apply(ctab,1,function(x){z=!is.na(x) & !is.nan(x); x[z] = x[z]/sum(x[z]); x[!z]=0; x}))
			wtab[is.element(wtab,0)] = NA
			jtab = apply(wtab*xtab,1,sum,na.rm=T)
			jtab[is.element(jtab,0)] = eval(parse(text=paste0("mean(u",xy,")")))
			XYwtd[[xy]] = jtab
		}
		## Rearrange data in format for `calcChloro'
		apple = apply(isumm,2,function(x){
			nday=length(x)
			slices = data.frame(X=XYwtd$X, Y=XYwtd$Y, julian=as.numeric(names(x)), zfld=x, region=rep(i,nday))
			names(slices) = gsub("zfld",zfld,names(slices))
			return(slices)
		})
		bushel = sapply(names(apple), function(x,a){
			bag = a[[x]]; nday=nrow(bag)
			bag = data.frame(slice=rep(x,nday),bag)
			return (bag)
		}, a=apple, simplify=FALSE) # TRUE craps out for some reason

		for(b in bushel)
			orchard = rbind(orchard,b)
		rownames(orchard) = 1:nrow(orchard)
	}
	names(orchard) = sub("slice",paste0(slice,collapse=delim),names(orchard))

	### Plot the data positions
	if (seePlot) {
		xlim = c(0,365)
		ylim = c(0,zmax)
		expandGraph(mfrow=c(1,1))
		plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="Julian Day",ylab="Chl a",xaxt="n")
		for (i in 1:length(crumble)) {
			ii = names(crumble)[i]
			iii = crumble[[i]]
			y  = attributes(iii)[[paste0(zfld,".mean")]];  x = as.numeric(names(y))
			lines(x,y,col=switch(i,"red","blue","green4"),lwd=2)
			points(x,y, pch=21, cex=1.2, col=switch(i,"red","blue","green4"), bg=switch(i,"pink","cyan","green"))
		}
		axis(1, at=x)
		legend("topright",lty=1,lwd=2,col=c("red","blue"),legend=anames,bty="n",inset=0.05)
	}
	### Save the data
	anam = paste0(onam,delim,"summary")
	omess = paste0(onam, "=orchard; ",
		"save(\"",onam,"\",file=\"",onam,".rda\"); ")
	omess = paste0(omess, anam,"=crumble; ",
		"save(\"",anam,"\",file=\"",anam,".rda\"); ")
	eval(parse(text=omess))
	return(invisible(orchard))
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~sliceDpolys


