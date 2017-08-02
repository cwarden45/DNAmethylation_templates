#setwd("R:\\Seq\\cwarden\\GENCODE_transcripts")

#uses filtered GTF created from filter_gtf.py script from https://github.com/cwarden45/RNAseq_templates/tree/master/Splicing_Workflow

gtf = "human/hg19/gencode.v19.annotation.FILTERED.gtf"
tss.flanking.length = 1500

bed = gsub(".gtf$",paste(".PROMOTER.",tss.flanking.length,"bp.bed",sep=""),gtf)

### ideally, edit above this line 						###
### warning: this will probably take a few hours to run ###


extract.id = function(char, id.type){
	char.info = unlist(strsplit(char,split=";"))
	info.text = char.info[grep(id.type, char.info)]
	info.text = gsub(id.type,"",info.text)
	info.text = gsub("\"","",info.text)
	info.text = gsub(" ","",info.text)
	return(info.text)
}#end def extract.id

gtf.table = read.table(gtf, head=F, sep="\t")
print(dim(gtf.table))
gtf.table = gtf.table[gtf.table$V3 == "transcript",]
print(dim(gtf.table))

transcriptID.values = sapply(as.character(gtf.table$V9), extract.id, id.type="transcript_id")
transcriptID.values = as.character(transcriptID.values)
transcriptID.values = gsub(".\\d+$","",transcriptID.values)
transcriptID = levels(as.factor(transcriptID.values))

geneName.values = sapply(as.character(gtf.table$V9), extract.id, id.type="gene_name")

gene.chr = gtf.table$V1[match(transcriptID, transcriptID.values)]
gene.strand = gtf.table$V7[match(transcriptID, transcriptID.values)]
gene.start = tapply(gtf.table$V4, transcriptID.values, min)
gene.stop = tapply(gtf.table$V4, transcriptID.values, max)

promoter.chr = gene.chr
promoter.start = c()
promoter.stop = c()
promoter.name = c()
promoter.score = c()
promoter.strand = gene.strand
for (i in 1:length(transcriptID)){
	tx.name = as.character(transcriptID[i])
	if(promoter.strand[i] == "+"){
		TSS = gene.start[i]
		promoter.start[i]=TSS-tss.flanking.length
		promoter.stop[i]=TSS+tss.flanking.length
	}else if(promoter.strand[i] == "-"){
		TSS = gene.stop[i]
		promoter.start[i]=TSS-tss.flanking.length
		promoter.stop[i]=TSS+tss.flanking.length
	}else {
		print(paste(promoter.strand[i]," is not a valid strand value",sep=""))
	}
	promoter.score[i]=0
	gene.name = as.character(geneName.values[match(tx.name,transcriptID.values)])
	promoter.name[i]=paste(gene.name,tx.name,"Promoter",sep="_")
}#end for (i in 1:nrow(txtable))

bed.table = data.frame(promoter.chr, promoter.start, promoter.stop,
						promoter.name, promoter.score, promoter.strand)
write.table(bed.table, bed, sep="\t", row.names=F, col.names=F, quote=F)