#===============================================================================
# Module 2: Samples
# -----------------
#  boxSeason.......Show seasonal data using boxplots.
#  plotDiversity...Plot diversity as barplots and overlay species richness.
#  trackComp.......Track the composition of phytoplankton groups
#===============================================================================


#boxSeason------------------------------2016-04-25
# Show seasonal data using boxplots.
#-----------------------------------------------RH
boxSeason <- function (fqtName="Ex01_Sample_Info", dbName="Examples.mdb",
     type="MDB", path=getwd(), fld="S", brks="M", pdf=FALSE, wmf=FALSE, ioenv=.GlobalEnv) {

	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))

	fnam=paste("Season-",fld,"x",brks,sep="")
	assign("PBSpton",list(module="P02_Temporal",call=match.call(),ioenv=ioenv,plotname=fnam),envir=.PBSptonEnv)
	if (type=="FILE") {
		eval(parse(text=paste("getFile(",fqtName,",senv=ioenv,try.all.frames=TRUE,tenv=penv()); dat=",fqtName,sep=""))) }
	else {
		getData(fqtName,dbName,type=type,path=path,tenv=penv()) 
		dat=PBSdat }
	ramp <- colorRamp(c("yellow","red","darkblue"))

	if (any(names(dat)=="YMD")) dat$date=as.POSIXct(dat$YMD)
	xrng <- range(dat$date)
	yr <- sort(unique(as.numeric(substring(dat$date,1,4))))
	yr = yr[1]:yr[length(yr)]; nyr <- length(yr)
	if (brks=="Q") { # Quarter
		moda0="01-01"; vlin=4
		moda=c("03-31","06-30","09-30","12-31"); 
		labs=c("Q1","Q2","Q3","Q4") }
	if (brks=="S") { # Season
		moda0="10-31"; vlin=4
		moda=c("02-28","05-31","08-31","10-31"); 
		labs=c("SP","SU","FA","WI") }
	if (brks=="B") { # Bimonthly
		moda0="01-01"; vlin=6
		moda=c("02-28","04-30","06-30","08-31","10-31","12-31"); 
		labs=c("JF","MA","MJ","JA","SO","ND") }
	if (brks=="M") { # Monthly
		moda0="01-01"; vlin=12
		moda=c("01-31","02-28","03-31","04-30","05-31","06-30","07-31","08-31","09-30","10-31","11-30","12-31"); 
		labs=c("Ja","Fe","Ma","Ap","My","Jn","Jl","Au","Se","Oc","No","De") }
	nmoda=length(moda)
	bcol=rgb(ramp(seq(0,1,len=nmoda)),maxColorValue=255)
	brks=paste(rep(yr,each=nmoda),rep(moda,nyr),sep="-")
	labs=paste(rep(yr,each=nmoda),rep(labs,nyr),sep="-")
	brks=c(paste(yr[1]-ifelse(moda0>moda[1],1,0),moda0,sep="-"),brks)

	dat$per <- cut(dat$date,as.POSIXct(brks),labels=labs)
	bag     <- split(dat[,fld],dat$per); nbag <- length(bag)
	xy      <- boxplot(bag,plot=FALSE)
	zna     <- apply(xy$stats,2,function(x){all(is.na(x))})
	ztrunc  <- clipVector(zna,"TRUE")
	ii      <- as.numeric(names(ztrunc)) # index of time periods which show data
#browser();return()
	ylim=c(0, max(sapply(bag[ii],function(x){if(length(x)>0) max(x) else 0 })))
	stuff=c("dat","brks","labs","bag","ii","bcol")
	packList(stuff,"PBSpton",tenv=.PBSptonEnv)

	if (pdf) pdf(file=paste(fnam,".pdf",sep=""),width=11,height=8.5)
	else if (wmf && .Platform$OS.type=="windows")
		do.call("win.metafile",list(filename=paste(fnam,".wmf",sep=""),width=11,height=8.5))
	else resetGraph()
	expandGraph(mar=c(6,5,2,1.5),mgp=c(2.75,.75,0))
	boxplot(bag[ii],las=2,range=0,lty=1,boxwex=0.5,staplewex=0,
		col=rep(bcol,length(ii))[ii],cex.lab=1.2,
		ylim=ylim, xlab=if(length(ii)==1) labs[!zna] else "",  # potentially problematic
		ylab=ifelse(fld=="S","Number of Species",ifelse(fld=="H","Shannon-Wiener Diversity Index","Evenness of Species")))
	zv=is.element(ii,seq(vlin,nbag,vlin))
	iv=(1:length(zv))[zv]+.5 # possible indices where vertical lines are drawn
	abline(v=iv,col="slategrey",lty=8)
	text(iv,par()$usr[4],yr[1:length(iv)],cex=1,adj=c(1.5,1.5))
	if (pdf|wmf) dev.off()
	invisible() }
#----------------------------------------boxSeason


#plotDiversity--------------------------2016-04-25
# Plot diversity as barplots and overlay species richness.
#-----------------------------------------------RH
plotDiversity <- function (fqtName="Ex01_Sample_Info",dbName="Examples.mdb",
     type="MDB",path=getwd(),bars="H", pnts="S", xnames=c("SID","Batch"),
     xint=40,clrs=c("skyblue","gold","blue","green4"), addlowess=TRUE, f=0.2, 
     pdf=FALSE, wmf=FALSE, ioenv=.GlobalEnv) {

	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))

	fnam=paste("Diversity-",bars,"x",pnts,sep="")
	assign("PBSpton",list(module="P02_Temporal",call=match.call(),ioenv=ioenv,plotname=fnam),envir=.PBSptonEnv)
	if (type=="FILE") {
		eval(parse(text=paste("getFile(",fqtName,",senv=ioenv,try.all.frames=TRUE,tenv=penv()); dat=",fqtName,sep=""))) }
	else {
		getData(fqtName,dbName,type=type,path=path,tenv=penv()) 
		dat=PBSdat }
	index=c("Shannon-Wiener diversity H", "Species richness S", "Species evenness E")
	names(index)=c("H","S","E")

	ypick=function(y) {
		ylim=extendrange(y,f=0.05); ylim[1]=0
		if (ylim[2]<1) fac=10^(-floor(log10(ylim[2]))) else fac=1
		ylim=ylim*fac
		ypos=ceiling(ylim[1]):floor(ylim[2]); ny=length(ypos)
		ypos=ypos[seq(1,ny,round(ny/min(6,ny)))]
		dfy=diff(ypos)[1]
		yint=diff(pretty(ypos[1]:dfy))[1]
		ytck=seq(ylim[1],ypos[length(ypos)],yint) #diff(ypos)[1]/4	)
		list(ylim=ylim/fac,ypos=ypos/fac,ytck=ytck/fac)
	}

	if (any(names(dat)=="YMD")) dat$date=as.POSIXct(dat$YMD)
	xrng <- range(dat$date)
	#yr <- sort(unique(as.numeric(substring(dat$YMD,1,4)))); nyr <- length(yr)
	xnams=apply(dat[,xnames,drop=FALSE],1,paste,collapse="-")
	nx=length(xnams)
	znams=seq(1,nx,round(nx/min(xint,nx)))
	xshow=rep("",nx); xshow[znams]=xnams[znams]; xchar=max(nchar(xshow))
	unpackList(ypick(dat[,bars]),scope="L")
	
	#if (eps) postscript(file=paste(fnam,".eps",sep=""),width=11,height=8.5,paper="special")
	if (pdf) pdf(file=paste(fnam,".pdf",sep=""),width=11,height=8.5)
	else if (wmf && .Platform$OS.type=="windows")
		do.call("win.metafile",list(filename=paste(fnam,".wmf",sep=""),width=11,height=8.5))
	else resetGraph()
	expandGraph(mar=c(max(3.5,ceiling(xchar^.7)),3.75,1,3.75),mgp=c(2.75,.75,0))
	x=barplot(height=dat[,bars],names.arg=xshow,space=0,col=clrs[1],cex.names=0.7,
		las=2,mgp=c(0,0.2,0),tck=.01,yaxt="n",xaxs="i",ylim=ylim)
	if (addlowess) lines(lowess(x,dat[,bars],f=f),col=clrs[3],lwd=2)
	axis(2,at=ytck,tck=-.005,labels=FALSE); axis(2,at=ypos,las=1)
	mtext(ifelse(any(bars==c("H","S","E")),index[bars],bars),side=2,line=2,cex=1.2)
	mtext("Sample identification",side=1,line=par()$mar[1]-1.5,cex=1.2)
	xlim=par()$usr[1:2]

	par(new=TRUE)
	unpackList(ypick(dat[,pnts]),scope="L")
	plot(x,dat[,pnts],type="n",xlim=xlim,xaxs="i",ylim=ylim,yaxs="i",axes=FALSE,xlab="",ylab="")
	abline(h=ypos[ypos>0],col="grey",lty=ifelse(pdf,1,3))
	if (addlowess) lines(lowess(x,dat[,pnts],f=f),col=clrs[4],lwd=2)
	points(x,dat[,pnts],pch=21,bg=clrs[2])
	axis(4,at=ytck,tck=-.005,labels=FALSE); axis(4,at=ypos,las=1)
#browser();return()
	mtext(ifelse(any(pnts==c("H","S","E")),index[pnts],pnts),side=4,line=1.75,cex=1.2)

	legend("topleft",pch=c(22,21),pt.bg=clrs[1:2],legend=c(bars,pnts),
		text.width=diff(par()$usr[1:2])*ifelse(addlowess,.10,.05),bg="aliceblue")
		if (addlowess) addLegend(.06,1,lty=c(1,1),col=clrs[3:4],legend=c("",""),bty="n")
	box()
	if (pdf|wmf) dev.off()
	packList(c("dat","x","xnams","xshow"),"PBSpton",tenv=.PBSptonEnv)
	invisible() }
#------------------------------------plotDiversity


#trackComp------------------------------2016-04-25
# Track the composition of phytoplankton groups
#-----------------------------------------------RH
trackComp = function(fqtName=c("Ex01_Sample_Info","Ex02_Species_Abundance"),
     dbName="Examples.mdb", type="MDB", path=getwd(), ndays=15, groups=NULL, 
     dlim=c("2007-01-01",substring(Sys.time(),1,19)),
     clrs=c("green","dodgerblue","orange","yellow","red","grey"),
     pdf=FALSE, wmf=FALSE, ioenv=.GlobalEnv, ...) {

	par0 <- par(no.readonly = TRUE)
	on.exit(par(par0))
	dots=list(...)
	fnam=paste("Comp",paste(dots$co,collapse="-"),"-Every-",ndays,"days",sep="")
	assign("PBSpton",list(module="P02_Temporal",call=match.call(),ioenv=ioenv,plotname=fnam),envir=.PBSptonEnv)
	if (type=="FILE") {
		eval(parse(text=paste("getFile(",fqtName[1],",senv=ioenv,try.all.frames=TRUE,tenv=penv()); dat=",fqtName[1],sep=""))) }
	else {
		eval(parse(text=paste("getData(\"",fqtName[1],"\",\"",dbName,"\",type=\"",type,"\",path=\"",path,"\",tenv=penv())",sep="")))
		dat=PBSdat }
	if (any(names(dat)=="YMD")) dat$date=as.POSIXct(dat$YMD)
	flds=names(dat)
	up=is.element(flds,c("PID","SID","EID","POS","YMD"))
	flds[!up]=tolower(flds[!up])
	names(dat)=flds
	dlim = as.POSIXct(dlim)
#browser();return()
	dat  = dat[dat$date>=dlim[1] & dat$date<=dlim[2] & !is.na(dat$date),]
	if (nrow(dat)==0) showError("Choose different date limits")
	if (any(flds=="dominant")) dat=dat[dat$dominant!="Zilcho" & dat$dominant!="" & !is.na(dat$dominant),]
	if (length(dots)>0) {
		unpackList(dots,scope="L")
		for (i in dots)
			eval(parse(text=paste("dat=biteData(dat,",i,")",sep="")))
	}
	yper=convYP(dat$date,ndays)
	dat$x=yper; dat$lab=names(yper)
	
	#Group proportions
	if (type=="FILE") {
		eval(parse(text=paste("getFile(",fqtName[2],",senv=ioenv,try.all.frames=TRUE,tenv=penv()); gdat=",fqtName[2],sep=""))) }
	else {
		eval(parse(text=paste("getData(\"",fqtName[2],"\",\"",dbName,"\",type=\"",type,"\",path=\"",path,"\",tenv=penv())",sep="")))
		gdat=PBSdat }
	gdat = gdat[is.element(gdat$SID,sort(unique(dat$SID))),]
	if (is.null(groups))
		groups=c("Skeletonema","Thalassiosira","Chaetoceros","Other Diatoms","Phytoflagellates","Grazers")
	px=sort(sapply(split(yper,names(yper)),unique))
	plob=list(x=px, lab=names(px)) # plot object
	for (i in groups) {
		idat=gdat[is.element(gdat$Group,i),]
		gA=idat$gA; names(gA)=idat$SID
		dat[,i] = gA[as.character(dat$SID)] 
		dat[is.na(dat[,i]),i]=0 
		plob[[i]]=sapply(split(dat[,i],dat$x),mean)  }
	tmp=as.data.frame(plob[groups])
	dfpic = data.frame(x=plob$x,base=rep(0,length(plob$x)),t(apply(tmp,1,cumsum)))
#browser();return()
	
	if (pdf) pdf(file=paste(fnam,".pdf",sep=""),width=11,height=8.5)
	else if (wmf && .Platform$OS.type=="windows")
		do.call("win.metafile",list(filename=paste(fnam,".wmf",sep=""),width=11,height=8.5))
	else resetGraph()
	expandGraph(mar=c(4.5,3.5,0.5,1.25),oma=c(0,0,0,0),las=1,fig=c(0,1,0,.95))
	plot(dfpic$x,dfpic$base,ylim=c(0,100),xlab="",ylab="",type="l",xaxs="i",yaxs="i",xaxt="n")
	dfgroups=gsub(" ",".",groups)
	for(i in 1:length(dfgroups)) {
		ii=dfgroups[i];  iii=ifelse(i==1,"base",dfgroups[i-1])
		xpol=c(dfpic$x,rev(dfpic$x))
		ypol=c(dfpic[,ii],rev(dfpic[,iii]))
		polygon(xpol,ypol,col=clrs[i],border=clrs[i])
	}
	abline(h=seq(10,90,10),v=sort(unique(round(dfpic$x))),col="moccasin",lty=3)
	axis(1,at=dfpic$x,labels=row.names(dfpic),las=2,tck=.01,mgp=c(0,.2,0),cex.axis=.8)
	mtext(paste("Year-period (interval =",ndays,"days)"),side=1,las=1,line=3,cex=1)
	mtext("Relative Composition (%)",side=2,las=0,line=2,cex=1.2)
	box(); 
	packList(c("dat","gdat","plob","dfpic"),"PBSpton",tenv=.PBSptonEnv)
	#Legend
	#par(new=TRUE,mar=c(0,2,1.25,0),fig=c(0,1,.95,1))
	par(new=TRUE,mar=c(0,1,0.1,0),fig=c(0,1,.95,1))
	plot(0,0,type="n",axes=FALSE,xlab="",ylab="")
	addLegend(.05,.5,fill=clrs,legend=groups,horiz=TRUE,bty="n",yjust=.5)
	if (pdf|wmf) dev.off()
	invisible() }
#----------------------------------------trackComp

