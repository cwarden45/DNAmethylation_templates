param.table = read.table("parameters.txt", header=T, sep="\t")
beta.prefix  = as.character(param.table$Value[param.table$Parameter == "beta_prefix"])
beta.normalization  = as.character(param.table$Value[param.table$Parameter == "beta_normalization"])

#downloaded from https://code.google.com/archive/p/bmiq/downloads
source("BMIQ_1.3.R")

#assumes you are working with EPIC data.  Need to modify code for other arrays.
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
data("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
annotation.table = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

probe.type = annotation.table$Type
probe.type[probe.type=="I"]=1
probe.type[probe.type=="II"]=2
names(probe.type)=annotation.table$Name

#Use parameters to specify beta table (parameters will then need to be changed after normalization)
output.normalization = paste(beta.normalization,"BMIQ",sep="_")

beta.in.file = paste(beta.prefix,"_",beta.normalization,".txt",sep="")

beta.table = read.table(beta.in.file, head=T, sep="\t")
SiteID = beta.table$SiteID
beta.mat = beta.table[,2:ncol(beta.table)]
rm(beta.table)

print(head(beta.mat))

for (i in 1:ncol(beta.mat)){
	temp.beta = as.numeric(beta.mat[,i])
	names(temp.beta) = SiteID
	
	temp.type = probe.type[!is.na(temp.beta)]
	temp.beta = temp.beta[!is.na(temp.beta)]
	
	bmiqObj = BMIQ(beta.v=temp.beta, design.v=temp.type,
					plots=TRUE, sampleID=names(beta.mat)[i])
	updated.beta = bmiqObj$nbeta
	updated.beta = as.numeric(updated.beta[match(SiteID, names(updated.beta))])
	updated.beta = round(updated.beta, digits=2)
	beta.mat[,i]=updated.beta
}#end for (i in ncol(beta.mat))

print(head(beta.mat))

output.table = data.frame(SiteID, beta.mat)
output.file = paste(beta.prefix,"_",output.normalization,".txt",sep="")
write.table(output.table, output.file, quote=F, sep="\t", row.names=F)