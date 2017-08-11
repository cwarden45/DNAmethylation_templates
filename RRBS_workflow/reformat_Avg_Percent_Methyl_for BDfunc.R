in.file = "TxDb_hg19_GENE_promoter_1500bp_average_methylation_Bismark_10x.txt"
out.file = "BDfunc/Avg_GENE_[cat].txt"
extract.gene.flag=FALSE

#in.file = "TxDb_hg19_CpG_ISLAND_IN_PROMOTER_1500bp_average_methylation_Bismark_10x.txt"
#out.file = "BDfunc/Avg_CpG_ISLAND_IN_PROMOTER_[cat].txt"
#extract.gene.flag=TRUE

meta.file = "sample_description.txt"
groupID = "[cat]"

extract.gene.from.region = function(char){
	char.info = unlist(strsplit(char,split="_"))
	return(char.info[1])
}#end def extract.gene.from.region

meta.table = read.table(meta.file,head=T,sep="\t")
meta.table = meta.table[!is.na(meta.table[,groupID]),]
sampleID = meta.table$userID
group = meta.table[,groupID]

percent.table = read.table(in.file,head=T,sep="\t")
percent.mat = round(percent.table[,match(sampleID, names(percent.table))], digits=2)

gene=percent.table[,1]
if(extract.gene.flag){
	gene=unlist(sapply(as.character(gene),extract.gene.from.region))
}#end if(extract.gene.flag)

output.table = data.frame(gene, percent.mat)
colnames(output.table) = c("gene",as.character(group))
write.table(output.table, out.file, sep="\t", row.names=F, quote=F)

#manually change header...
#BD-Func
#groupA	groupB groupA groupB groupA groupB