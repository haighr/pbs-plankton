make.phytospp = function()
{
	getData("phytospp", "PBSplankton", type="MDB")
	phytospp = PBSdat
	colnames(phytospp) = tolower(colnames(phytospp))
	phytospp = phytospp[order(phytospp$pid,phytospp$sid),]

	getData("C_Toxins", "PBSplankton", type="MDB")
	toxins = PBSdat
	colnames(toxins) = tolower(colnames(toxins))
	toxins = toxins[order(toxins$acronym),]
	eff = (match("toxin",colnames(toxins))+1):ncol(toxins)  ## index of toxin effects
	toxins[,eff][is.na(toxins[,eff])] = 0
	toxins[,eff] = sapply(toxins[,eff],as.logical)
	attr(phytospp, "toxins") = toxins

	stamp = Sys.time()
	#today = gsub("-","",substring(stamp,3,10))
	attr(phytospp,"created") = stamp
	if (file.exists("phytospp.rda")) {
		old.stamp = file.info("phytospp.rda")[,"mtime"]
		old.date  = gsub("-","",substring(old.stamp,3,10))
		file.copy(from="phytospp.rda", to=paste0("phytospp-",old.date,".rda"), overwrite=TRUE, copy.date=TRUE)
	}
	save("phytospp", file="phytospp.rda")
#browser();return()
}

make.phytospp()
