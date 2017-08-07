param.table = read.table("parameters_methylKit.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
min.cov=as.numeric(as.character(param.table$Value[param.table$Parameter == "Min_Coverage"]))
quant.type=as.character(param.table$Value[param.table$Parameter == "Quantification_Method"])
output.prefix=as.character(param.table$Value[param.table$Parameter == "methyl_percent_prefix"])
dmr.vars  = as.character(param.table$Value[param.table$Parameter == "dmr_groups"])
output.format  = as.character(param.table$Value[param.table$Parameter == "COHCAP_output_format"])
overall.result.folder = as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
project.folder = as.character(param.table$Value[param.table$Parameter == "project_folder"])
project.name = as.character(param.table$Value[param.table$Parameter == "comp_name"])
expression.file = as.character(param.table$Value[param.table$Parameter == "expression_file"])
promoter.bed = as.character(param.table$Value[param.table$Parameter == "gene_region_BED"])
COHCAP.mapping = as.character(param.table$Value[param.table$Parameter == "gene_region_mapping"])
wig.output = as.character(param.table$Value[param.table$Parameter == "wig_types"])
num.groups = as.character(param.table$Value[param.table$Parameter == "COHCAP_num_groups"])
pair.var = as.character(param.table$Value[param.table$Parameter == "COHCAP_paired"])
site.pvalue = as.numeric(as.character(param.table$Value[param.table$Parameter == "site_pvalue"]))
site.fdr = as.numeric(as.character(param.table$Value[param.table$Parameter == "site_fdr"]))
site.delta.beta = as.numeric(as.character(param.table$Value[param.table$Parameter == "site_delta_beta"]))
island.pvalue = as.numeric(as.character(param.table$Value[param.table$Parameter == "island_pvalue"]))
island.fdr = as.numeric(as.character(param.table$Value[param.table$Parameter == "island_fdr"]))
island.delta.beta = as.numeric(as.character(param.table$Value[param.table$Parameter == "island_delta_beta"]))
max.cluster.dist = as.numeric(as.character(param.table$Value[param.table$Parameter == "max_cluster_dist"]))
sites.per.island = as.numeric(as.character(param.table$Value[param.table$Parameter == "min_sites_per_island"]))
methyl.threshold = as.numeric(as.character(param.table$Value[param.table$Parameter == "methyl_cutoff"]))
unmethyl.threshold = as.numeric(as.character(param.table$Value[param.table$Parameter == "unmethyl_cutoff"]))
trt.group = as.character(param.table$Value[param.table$Parameter == "treatment_group"])
library.type = as.character(param.table$Value[param.table$Parameter == "Read_Pairing"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])

sep.cov.folder = paste(overall.result.folder,"/Bismark_Coverage_Files",sep="")
project.folder = paste(overall.result.folder,"/",project.folder,sep="")
dir.create(project.folder)
dmr.file = paste(project.folder,"/",project.name,"_methylKit_DMR_stats.txt",sep="")
region.percent.file = paste(project.folder,"/",project.name,"_methylKit_region_percent_methylated.txt",sep="")

library(methylKit)

bismark.cov.files = list.files(sep.cov.folder)
bismark.cov.files = list.files(sep.cov.folder)
bismark.cov.files = bismark.cov.files[grep(".bismark.cov",bismark.cov.files)]

if(library.type == "PE"){
	sampleID = gsub("_S\\d+_L\\d{3}_R1_001_val_1_bismark_bt2_pe.bismark.cov","",bismark.cov.files)
}else{
	stop("Add regular expression for Bismark SE reads")
}

bismark.cov.files = paste(sep.cov.folder,"/",bismark.cov.files,sep="")


sample.table = read.table(sample.description.file, sep="\t", header=T)
userID = as.character(sample.table$userID)
dmr.vars = unlist(strsplit(dmr.vars,","))
dmr.cols = sample.table[,dmr.vars]

prim.2group = c()

if(is.null(dim(dmr.cols))){
	names(dmr.cols) = userID
	dmr.cols = dmr.cols[!is.na(dmr.cols)]
	userID = names(dmr.cols)
	if(length(levels(dmr.cols)) == 2){
		prim.2group = rep(0,length(dmr.cols))
		prim.2group[dmr.cols == trt.group]=1
		ref.group = levels(dmr.cols)[levels(dmr.cols) != trt.group]
	}else{
		stop("methylKit requires 2-group primary variable (with or without covariates)")
	}

}else{
	stop("Add code for multi-variate methylKit analysis")
	rownames(dmr.cols) = userID
	dmr.cols = na.omit(dmr.cols)
	userID = rownames(dmr.cols)
}

sample.table = sample.table[match(userID, as.character(sample.table$userID)),]
comp.files = bismark.cov.files[match(sample.table$sampleID,sampleID)]
methylKit.files =split(comp.files, 1:length(comp.files))
methylKit.ids = split(as.character(sample.table$userID), 1:length(sample.table$userID))
myobj = methRead(location=methylKit.files, sample.id=methylKit.ids, assembly=genome,
					pipeline="bismarkCoverage", mincov=min.cov, dbdir=project.folder,
					treatment=prim.2group)
meth=unite(myobj)

promoter.table = read.table(promoter.bed, head=F, sep="\t")
print(dim(promoter.table))
if(length(grep("_",promoter.table$V1)) > 0){
	promoter.table = promoter.table[-grep("_",promoter.table$V1),]
	promoter.table$V1 = as.factor(as.character(promoter.table$V1))
	print(dim(promoter.table))
}
refGR = GRanges(Rle(promoter.table$V1),
				IRanges(start=promoter.table$V2, end=promoter.table$V3),
				Names=promoter.table$V4,
				Rle(strand(promoter.table$V6)))
refGR=unique(refGR)
GR.table = data.frame(refGR)
GR.ID = paste(GR.table$seqnames,":",GR.table$start,"-",GR.table$end,sep="")
					
regions=regionCounts(meth, regions=refGR)
region.table = data.frame(regions)
percent.regionID = paste(region.table$chr,":",region.table$start,"-",region.table$end,sep="")
for(i in 0:(length(sample.table$userID)-1)){
	cov.index = 5 + i * 3
	meth.index = 5 + i * 3 + 1
	total.counts = as.numeric(region.table[,cov.index])
	methyl.counts = as.numeric(region.table[,meth.index])
	sample.percent = 100 * methyl.counts / total.counts
	
	if(i ==0){
		percent.table = data.frame(sample.percent)
	}else{
		percent.table = data.frame(percent.table, sample.percent)
	}
}#end for(i in 1:length(sample.table$userID))
names(percent.table)=sample.table$userID

avgGroupExpression = function (geneExpr, groups) {
	avg.expr = tapply(geneExpr, groups, mean, na.rm=T)
	return(avg.expr)
}#end def avgGroupExpression

avg.percent.methyl = t(apply(percent.table, 1, avgGroupExpression, groups=prim.2group))
colnames(avg.percent.methyl)=paste(c(ref.group, trt.group),".avg.percent.methylation",sep="")
regionID = GR.table$Names[match(percent.regionID, GR.ID)]
region.percent.table = data.frame(regionID, region.table[,1:4], percent.table)
write.table(region.percent.table, region.percent.file, quote=F, sep="\t", row.names=F)


myDiff=data.frame(calculateDiffMeth(regions))
status = rep("No Change",nrow(myDiff))
status[(myDiff$meth.diff > island.delta.beta)&(myDiff$pvalue < island.pvalue)&(myDiff$qvalue < island.fdr)]=paste(trt.group," Increased Methylation",sep="")
status[(myDiff$meth.diff < -island.delta.beta)&(myDiff$pvalue < island.pvalue)&(myDiff$qvalue < island.fdr)]=paste(trt.group," Decreased Methylation",sep="")
print(table(status))

myDiff.ID = paste(myDiff$chr,":",myDiff$start,"-",myDiff$end,sep="")

#if different transcripts have same TSS, there can be fewer names than in the original table, but gene name should be the same
regionID = GR.table$Names[match(myDiff.ID, GR.ID)]

avg.percent.methyl = avg.percent.methyl[match(myDiff.ID,percent.regionID),]

output.table = data.frame(regionID, avg.percent.methyl, myDiff, status)
write.table(output.table, dmr.file, quote=F, sep="\t", row.names=F)
