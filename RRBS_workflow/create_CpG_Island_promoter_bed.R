island.name = "UCSC_CpG_Islands"
CpG.Island.bed = "UCSC_CpG_Islands.bed"
genome = "hg19"
tss.flanking.length = 1500

library("IRanges")

promoter.bed = paste("TxDb_",genome,"_promoter_",tss.flanking.length,"bp.bed",sep="")
output.bed = paste("TxDb_",genome,"_",island.name,"_in_promoter_",tss.flanking.length,"bp.bed",sep="")

overlap.bed = "temp_overlap.bed"

command = paste("/opt/bedtools2/bin/bedtools intersect -wb -a ",CpG.Island.bed,
				" -b ",promoter.bed, " > ",overlap.bed,sep="")
system(command)

overlap.table = read.table(overlap.bed, head=F, sep="\t")

extract.gene.region = function(char){
	char.info = unlist(strsplit(char,split="_"))
	gene = char.info[1]
	return(gene)
}#end def extract.gene.region
	
gene.symbol = as.character(sapply(as.character(overlap.table$V8),extract.gene.region))
overlap.table$V4 = as.character(overlap.table$V4)
overlap.table$V4 = gsub(":_","",overlap.table$V4)
newName = paste(gene.symbol,overlap.table$V4,sep="_")

overlap.table$V2 = as.numeric(overlap.table$V2)
overlap.table$V3 = as.numeric(overlap.table$V3)
output.table = data.frame(overlap.table$V1, overlap.table$V2, overlap.table$V3,
							newName, score=rep(0,length(newName)),overlap.table$V10)
print(dim(output.table))
table.txt = apply(output.table, 1, paste, collapse="\t")
reduce.regions=function(text.arr){
	mat = matrix(unlist(strsplit(text.arr,split="\t")),ncol=6, byrow=TRUE)
	region.chr = mat[,1]
	non.canonical = region.chr[grep("_",region.chr)]
	if ((length(non.canonical) > 0)&(length(non.canonical) != nrow(mat))){
		if(nrow(mat) - length(non.canonical) == 1){
			mat = mat[-grep("_",region.chr),]
			return.txt = paste(mat, collapse="\t")
			return(return.txt)
		}else{
			mat = mat[-grep("_",region.chr),]
			#print(mat)
			region.chr = mat[,1]
		}
	}else if (length(non.canonical) > 0){
		return(NA)
	}
	region.start = as.numeric(mat[,2])
	region.stop = as.numeric(mat[,3])
	
	ir = reduce(IRanges(start=region.start, end=region.stop))
	reduced.start = start(ir)
	reduced.stop = end(ir)
	all.pos = c(reduced.start,reduced.stop)
	
	return.txt = paste(mat[1,1],min(all.pos),max(all.pos),
						mat[1,4], 0, mat[1,6], sep="\t")
	return(return.txt)
}
output.txt = tapply(table.txt, newName, reduce.regions)
output.txt = as.character(output.txt[!is.na(output.txt)])
output.table = t(data.frame(strsplit(output.txt, split="\t")))
print(dim(output.table))
output.table = output.table[!is.na(output.table[,2]),]
print(dim(output.table))
write.table(output.table, output.bed, col.names=F, row.names=F, sep="\t", quote=F)

file.remove(overlap.bed)
