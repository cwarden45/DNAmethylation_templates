outputPrefix = "methylKit_heatmap"
plot.groups = c("Group")
cluster.distance = "Pearson_Dissimilarity"

min.percent.observed = 0.75
#methylKit should only consider fully covered regions, but keep this parameter to potentially clean up regions

library("gplots")

param.table = read.table("parameters_methylKit_CpG_Island.txt", header=T, sep="\t")
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

project.folder = paste(overall.result.folder,"/",project.folder,sep="")
heatmap.file = paste(outputPrefix,"_",project.name,".png",sep="")

dmr.file = paste(project.folder,"/",project.name,"_methylKit_DMR_stats.txt",sep="")
region.percent.file = paste(project.folder,"/",project.name,"_methylKit_region_percent_methylated.txt",sep="")

#reformat percent methylation
percent.table = read.table(region.percent.file, head=T, sep="\t")
#assumes region starts with gene symbol
extract.gene.from.region = function(char){
	char.info = unlist(strsplit(char,split="_"))
	return(char.info[1])
}#end def extract.gene.from.region

gene=as.character(unlist(sapply(as.character(percent.table$regionID), extract.gene.from.region)))
percent.mat = percent.table[,6:ncol(percent.table)]
percent.table = data.frame(island=percent.table$regionID, gene = gene,percent.mat)
print(dim(percent.table))
count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
}#end count.covered
covered.sample.count = apply(percent.mat, 1, count.covered)
percent.table = percent.table[covered.sample.count >= min.percent.observed * ncol(percent.mat),]
print(dim(percent.table))

dmr.table = read.table(dmr.file, head=T, sep="\t")
dmrs = dmr.table$regionID[dmr.table$status != "No Change"]

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

sample.description.table = read.table(sample.description.file, sep="\t", head=T)
sample.label = sample.description.table$userID

temp.percent = percent.table[match(dmrs, percent.table$island, nomatch=0),]
gene.labels = percent.table$gene[match(dmrs, percent.table$island, nomatch=0)]
temp.percent = temp.percent[,match(sample.label, names(percent.table))]
print(dim(temp.percent))

	cor.dist = function(mat){
		cor.mat = cor(as.matrix(t(mat)))
		dis.mat = 1 - cor.mat
		return(as.dist(dis.mat))	
	}#end def cor.dist

	if (cluster.distance == "Pearson_Dissimilarity"){
		print("Using Pearson Dissimilarity as Distance in Heatmap...")
		dist.fun = cor.dist
	}else{
		dist.fun=dist
	}

	if(length(plot.groups) > 1){
		source("heatmap.3.R")
		grp1 = as.character(sample.description.table[,plot.groups[1]])
		grp2 = as.character(sample.description.table[,plot.groups[2]])
		group.levels = c(levels(as.factor(grp1)),levels(as.factor(grp2)))

		color.palette <- fixed.color.palatte[1:length(group.levels)]
		labelColors1 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		labelColors2 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		
		ab.percent = t(temp.percent)
		if(length(gene.labels) < 25){
			col.labels = gene.labels
		} else {
			col.labels = rep("", length(gene.labels))
		}
		colnames(ab.percent) = rep("", length(gene.labels))
		rownames(ab.percent) = sample.label

		column_annotation <- as.matrix(gene.labels)
		colnames(column_annotation) <- c("")

		row_annotation <- data.frame(label1 = labelColors1, label2 = labelColors2)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) <- c(plot.groups)

		png(file = heatmap.file)
		heatmap.3(ab.percent,   distfun = dist.fun, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					RowSideColors=row_annotation, trace="none",labCol=col.labels,
					margins = c(8,13),RowSideColorsSize=4, dendrogram="both")
		legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		dev.off()
	} else {
		grp = sample.description.table[,plot.groups]
		group.levels = levels(as.factor(grp))

		color.palette <- fixed.color.palatte[1:length(group.levels)]
		labelColors = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors[grp == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))

		ab.percent = t(temp.percent)
		if(length(gene.labels) < 25){
			col.labels = gene.labels
		} else {
			col.labels = rep("", length(gene.labels))
		}
		colnames(ab.percent) = rep("", length(gene.labels))
		rownames(ab.percent) = sample.label
		
		png(file = heatmap.file)
		heatmap.2(ab.percent, distfun = dist.fun, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 RowSideColors=labelColors, trace="none", margins = c(5,15),labCol=col.labels)
		legend("topright", legend=group.levels,	col=color.palette, pch=15, cex=0.7)
		dev.off()
	}#end else