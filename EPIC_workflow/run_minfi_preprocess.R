param.table = read.table("parameters.txt", header=T, sep="\t")
description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
beta.prefix  = as.character(param.table$Value[param.table$Parameter == "beta_prefix"])
beta.normalization  = as.character(param.table$Value[param.table$Parameter == "beta_normalization"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])

library("minfi")

sample.table = read.table(sample.description.file, sep="\t", header=T)
sampleID = paste(sample.table$Sentrix_ID, sample.table$Sentrix_Position, sep="_")
userID = as.character(sample.table$userID)

meta.table = read.table(description.file, head=T, sep="\t")
RG.raw = read.metharray(basenames=meta.table$Base)

if (beta.normalization == "illumina"){
	methyl.norm = preprocessIllumina(RG.raw, bg.correct = TRUE, normalize = "controls", reference = 1)
	beta.table = getBeta(methyl.norm)
	beta.table = beta.table[,match(sampleID, colnames(beta.table))]
	probes = rownames(beta.table)
	colnames(beta.table) = userID
	
	output.table = data.frame(SiteID=probes, beta.table)
	beta.file = paste(beta.prefix,"_",beta.normalization,".txt",sep="")
	write.table(output.table, file=beta.file, sep="\t", quote=F, row.names=F)
} else if (beta.normalization == "funnorm"){
	methyl.norm = preprocessFunnorm(RG.raw, nPCs=2, sex = NULL, bgCorr = TRUE, dyeCorr = TRUE)
	beta.table = getBeta(methyl.norm)
	beta.table = beta.table[,match(sampleID, colnames(beta.table))]
	probes = rownames(beta.table)
	colnames(beta.table) = userID

	output.table = data.frame(SiteID=probes, beta.table)
	beta.file = paste(beta.prefix,"_",beta.normalization,".txt",sep="")
	write.table(output.table, file=beta.file, sep="\t", quote=F, row.names=F)
} else if (beta.normalization == "noob"){
	methyl.norm = preprocessNoob(RG.raw, offset = 15, dyeCorr = TRUE)
	beta.table = getBeta(methyl.norm)
	beta.table = beta.table[,match(sampleID, colnames(beta.table))]
	probes = rownames(beta.table)
	colnames(beta.table) = userID

	output.table = data.frame(SiteID=probes, beta.table)
	beta.file = paste(beta.prefix,"_",beta.normalization,".txt",sep="")
	write.table(output.table, file=beta.file, sep="\t", quote=F, row.names=F)
}else{
	stop("beta_normalization must be \"illumina\", \"funnorm\", or \"noob\"")
}