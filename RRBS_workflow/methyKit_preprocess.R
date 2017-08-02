param.table = read.table("parameters.txt", header=T, sep="\t")
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
min.cov=as.numeric(param.table$Value[param.table$Parameter == "Min_Coverage"])
alignment.folder=as.character(param.table$Value[param.table$Parameter == "Alignment_Folder"])
library.type=as.character(param.table$Value[param.table$Parameter == "Read_Pairing"])

check.overlap=FALSE
if(library.type == "PE"){
	#probably better to use Bismark quantification for PE data, but you can compare
	check.overlap=TRUE
}#end if(library.type == "PE")

library("methylKit")

bam.files = list.files(alignment.folder)
bam.files = bam.files[grep(".bam$",bam.files)]

sort.prefiles = bam.files[grep("\\.\\d{4}\\.bam",bam.files)]
if(length(sort.prefiles) > 0){
	bam.files = bam.files[-grep("\\.\\d{4}\\.bam",bam.files)]
}

for(bam.file in bam.files){
	sampleID = gsub(".bam","",bam.file)
	count.table = paste(alignment.folder,"/",sampleID,"/",sampleID,"_CpG.txt",sep="")
	if(!(file.exists(count.table))){
		print(paste("Quantifying counts for ",sampleID,"...",sep=""))
		output.folder = paste(alignment.folder,"/",sampleID,sep="")
		fullPath = paste(alignment.folder,"/",bam.file,sep="")
		methRaw = processBismarkAln(location=fullPath, sample.id=sampleID,
									nolap=check.overlap,save.context="CpG",
									mincov=min.cov, assembly=genome, save.folder=output.folder)
	}else{
		print(paste("File already created for ",sampleID,"...",sep=""))
	}
}#end for(bam.file in bam.files)