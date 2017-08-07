#assumes you set COHCAP to output files in ".txt" format

compID = ""
overall.result.folder = "../Result"
project.folder = "COHCAP_Results"
metadata.table = "sample_description.txt"
min.percent.observed = 0.75
cluster.distance="Pearson Dissimilarity"
plot.groups = c("Group")

heatmap.file  = paste(compID,"_island_heatmap_v2.png",sep="")
project.folder = paste(overall.result.folder,"/",project.folder,sep="")
island.beta.values = paste(project.folder,"/Raw_Data/Promoter_",compID,"_CpG_island_filtered_beta_values-Avg_by_Island.txt",sep="")
island.dmr.list = paste(project.folder,"/CpG_Island/Promoter_",compID,"_CpG_island_filtered-Avg_by_Island.txt",sep="")


###edit code above this line###

library(gplots)

dmr.table = read.table(island.dmr.list, sep="\t", head=T)
dmr.regions = dmr.table$island
dmr.genes = dmr.table$gene
	
fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black","darkgreen","thistle1","tan","orchid1",colors())

sample.description.table = read.table(metadata.table, sep="\t", head=T)
sample.label = sample.description.table$userID
	
beta.table = read.table(island.beta.values, sep="\t", head=T)
beta.islands = as.character(beta.table$island)
	
temp.beta = beta.table[match(dmr.regions, beta.islands, nomatch=0),]
rownames(temp.beta)=temp.beta$island
temp.beta = temp.beta[,match(sample.label, names(temp.beta))]

print(dim(temp.beta))
count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
}#end count.covered
covered.sample.count = apply(temp.beta, 1, count.covered)
temp.beta = temp.beta[covered.sample.count >= min.percent.observed * ncol(temp.beta),]
dmr.genes = dmr.genes[match(rownames(temp.beta),dmr.regions)]
print(dim(temp.beta))

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

		color.palette = fixed.color.palatte[1:length(group.levels)]
		labelColors1 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		labelColors2 = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))
		
		ab.methyl = t(temp.beta)
		if(length(dmr.genes) < 30){
			colnames(ab.methyl) = dmr.genes
		} else {
			colnames(ab.methyl) = rep("", length(dmr.genes))
		}
		rownames(ab.methyl) = sample.label

		column_annotation = as.matrix(dmr.genes)
		colnames(column_annotation) = c("")

		row_annotation = data.frame(label1 = labelColors1, label2 = labelColors2)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) = c(plot.groups)

		png(file = heatmap.file)
		heatmap.3(ab.methyl,   distfun = dist.fun, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					RowSideColors=row_annotation, trace="none", margins = c(8,13),RowSideColorsSize=4, dendrogram="both")
		legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		dev.off()
	} else {
		grp = sample.description.table[,plot.groups]
		group.levels = levels(as.factor(grp))

		color.palette = fixed.color.palatte[1:length(group.levels)]
		labelColors = rep("black",times=length(sample.label))
		for (i in 1:length(group.levels)){
			labelColors[grp == as.character(group.levels[i])] = color.palette[i]
		}#end for (i in 1:length(group.levels))

		ab.methyl = t(temp.beta)
		if(length(dmr.genes) < 30){
			colnames(ab.methyl) = dmr.genes
		} else {
			colnames(ab.methyl) = rep("", length(dmr.genes))
		}
		rownames(ab.methyl) = sample.label
		
		png(file = heatmap.file)
		heatmap.2(ab.methyl,   distfun = dist.fun, col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
					 RowSideColors=labelColors, trace="none", margins = c(5,15))
		legend("topright", legend=group.levels,	col=color.palette, pch=15, cex=0.7)
		dev.off()
	}#end else