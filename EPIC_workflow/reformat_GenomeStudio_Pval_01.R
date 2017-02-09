#Export Group beta values, defining each array as a group
FinalReport = "GenomeStudio_Files/FinalReport.txt"

parameterFile = "parameters_GenomeStudio.txt"

###### Should be able to just edit text above this line ######

#FYI, COHCAP has a couple default 450k island mappings (UCSC Islands or HMM Islands).
#However, you may still want to provide a custom mapping for TSS1500 + TSS200 + 5' UTR + 1stExon.
#This type of mapping file is available for 450k and EPIC array from https://sourceforge.net/projects/cohcap/files/additional_Bioconductor_annotations.zip

#The region-mapping doesn't matter for this script,
#	but it is something that you might want to consider when running COHCAP

param.table = read.table(parameterFile, header=T, sep="\t")
beta.prefix  = as.character(param.table$Value[param.table$Parameter == "beta_prefix"])
beta.normalization  = as.character(param.table$Value[param.table$Parameter == "beta_normalization"])
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])

input.table = read.table(FinalReport, head=T, sep="\t", skip=8)
print(dim(input.table))
probeID = input.table[,1]
#remove last, empty column
input.table = input.table[,-ncol(input.table)]

beta.mat = round(input.table[,seq(2,ncol(input.table),2)], digits=2)
pvalue.mat = input.table[,seq(3,ncol(input.table),2)]
rm(input.table)

filtered.beta.mat = beta.mat
filtered.beta.mat[pvalue.mat > 0.01] = NA
rm(beta.mat)
rm(pvalue.mat)

meta.table = read.table(sample.description.file,head=T, sep="\t")

sample.label = meta.table$userID[match(colnames(filtered.beta.mat),paste("X",meta.table$Sentrix_ID,"_",meta.table$Sentrix_Position,".AVG_Beta",sep=""))]
colnames(filtered.beta.mat) = sample.label

print("Writing COHCAP beta file to text")
beta.file = paste(beta.prefix,"_",beta.normalization,".txt",sep="")
cohcap.input = data.frame(SiteID = probeID, filtered.beta.mat)
write.table(cohcap.input, beta.file, sep="\t", row.names=F, quote=F)
