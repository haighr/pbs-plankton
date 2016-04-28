#===============================================================================
# Module 3: Satellites
# --------------------
#  calcChloro......Calculate chlorophyll peaks and integrate.
#  calcDeriv.......Calculate derivatives of loess-smoothe curve.
#===============================================================================


#calcChloro-----------------------------2016-04-28
# Calculate chlorophyll peaks and integrate
#-----------------------------------------------RH
calcChloro  = function(dat, wd=getwd(), xfld="period", yfld="Chlo", 
  xint=0.1, span=0.75, diff.required=FALSE, slice=NULL, days.per=8,
  salmEnd="Jun-15",
  show=c(fit=TRUE,der1=TRUE,der2=FALSE),
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
	satellite = strsplit(infile,split="\\.")[[1]][1]  ## assumes leading string before first "." is satellite name
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
			salmon = sapply(zbloo,function(x){matrix(x[,1],ncol=1)},simplify=FALSE)
			salmon = sapply(zbloo,function(x){data.frame(salmon=x[,1])},simplify=FALSE)
			salmon$x[2,1] = salmX1
			salmon$y[2,1] = predFunc(salmX1)
			salmBloo = sapply(salmon$x,function(x){stats::integrate(predFunc,x[1],x[2])$value})
			salmDur = salmX1 - zbloo$x[mbloo][1,1]
			#salmon = sapply(salmon,function(x){dimnames(x)=list(c("start","end"),"salmon");x},simplify=FALSE)
			#salmBloo = stats::integrate(predFunc,salmon$x[1],salmon$x[2])$value
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
			#abline(v=167,col="lightgrey")
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
	return(res.slice)
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

