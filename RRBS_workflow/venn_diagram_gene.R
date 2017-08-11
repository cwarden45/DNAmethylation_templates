library("Vennerable")
#install manually

cohcap.file = "../../Result/COHCAP_Results/CpG_Island/Promoter_[compID]_CpG_island_filtered-Avg_by_Island.txt"
methylkit.file = "../../Result/methylKit_Results/[compID]_CpG_Island_methylKit_DMR_stats.txt"

extract.gene.from.region = function(char){
	char.info = unlist(strsplit(char,split="_"))
	return(char.info[1])
}#end def extract.gene.from.region

cohcap.table = read.table(cohcap.file, head=T, sep="\t")
methylkit.table = read.table(methylkit.file, head=T, sep="\t")

cohcap.up = as.character(cohcap.table$gene[cohcap.table[,5] > 0])
cohcap.down = as.character(cohcap.table$gene[cohcap.table[,5] < 0])

methylkit.up = as.character(methylkit.table$regionID[methylkit.table$status == "Trt2 Increased Methylation"])
methylkit.up = unlist(sapply(as.character(methylkit.up),extract.gene.from.region))
methylkit.down = as.character(methylkit.table$regionID[methylkit.table$status == "Trt2 Decreased Methylation"])
methylkit.down = unlist(sapply(as.character(methylkit.down),extract.gene.from.region))

gene.list = list(C.UP=cohcap.up,
		C.DOWN=cohcap.down,
		M.UP=methylkit.up,
		M.DOWN=methylkit.down)
vennObj = Venn(gene.list)

png("DMR_venn.png")
#if 2 or 3 groups, use , doWeights = FALSE to remove scaling
plot(vennObj, type="ellipses")
dev.off()
