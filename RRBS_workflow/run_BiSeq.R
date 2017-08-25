min.percent.observed = 0.75

param.table = read.table("parameters_BiSeq.txt", header=T, sep="\t")
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
threads=as.numeric(as.character(param.table$Value[param.table$Parameter == "Threads"]))

sep.cov.folder = paste(overall.result.folder,"/Bismark_Coverage_Files",sep="")
sep.cov.folder = "../../Result/Bismark_Coverage_Files"
project.folder = paste(overall.result.folder,"/",project.folder,sep="")
dir.create(project.folder)
dms.file = paste(project.folder,"/",project.name,"_BiSeq_DMS_stats.txt",sep="")
dmr.file = paste(project.folder,"/",project.name,"_BiSeq_DMR_stats.txt",sep="")
site.percent.file = paste(project.folder,"/",project.name,"_BiSeq_site_percent_methylated.txt",sep="")
region.percent.file = paste(project.folder,"/",project.name,"_BiSeq_region_percent_methylated.txt",sep="")
RData.file = paste(project.folder,"/",project.name,"_BiSeq.RData",sep="")

library(BiSeq)

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

if(is.null(dim(dmr.cols))){
	names(dmr.cols) = userID
	dmr.cols = dmr.cols[!is.na(dmr.cols)]
	userID = names(dmr.cols)
	var1=dmr.cols
}else{
	stop("Add code for multi-variate BiSeq analysis")
	rownames(dmr.cols) = userID
	dmr.cols = na.omit(dmr.cols)
	userID = rownames(dmr.cols)
	var1=dmr.cols[,1]
}

if(trt.group != "continuous"){
	grp.levels = levels(as.factor(var1))
	#delta-beta is group1 - group2 (so, you want group1 to be the treatment)
	var1 = factor(var1, levels = c(trt.group,grp.levels[grp.levels != trt.group]))
}else{
	var1=as.numeric(var1)
}

sample.table = sample.table[match(userID, as.character(sample.table$userID)),]
comp.files = bismark.cov.files[match(sample.table$sampleID,sampleID)]

BSraw.obj = readBismark(comp.files, as.character(sample.table$sampleID))

print("Site Clustering")
if(trt.group != "continuous"){
	BSraw.obj= clusterSites(object = BSraw.obj,
						groups = var1, perc.samples = min.percent.observed,
						min.sites = sites.per.island, max.dist = max.cluster.dist)
}else{
	BSraw.obj= clusterSites(object = BSraw.obj,
						perc.samples = min.percent.observed,
						min.sites = sites.per.island, max.dist = max.cluster.dist)

}
print(BSraw.obj)
						
print("Import Region Information")
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

#print(head(promoter.table))
#print(refGR)

BSraw.obj =  subsetByOverlaps(BSraw.obj, refGR)
ov = findOverlaps(BSraw.obj, refGR)
rowRanges(BSraw.obj)$cluster.id[queryHits(ov)] = as.character(refGR$Names[subjectHits(ov)])
print(BSraw.obj)					
					
ind.cov = totalReads(BSraw.obj) > 0
quant = quantile(totalReads(BSraw.obj)[ind.cov], 0.9)

print(paste("Set maximum coverage of ",quant,"x",sep=""))
BSraw.obj = limitCov(BSraw.obj, maxCov = quant)
print(BSraw.obj)

print("Estimate Methylation Values")
predictedMeth = predictMeth(object = BSraw.obj, mc.cores=threads)

site.methyl = methLevel(predictedMeth)
clusterID = rowData(predictedMeth)
clusterID = clusterID$cluster.id
site.methyl = data.frame(clusterID, site.methyl)
write.table(site.methyl, site.percent.file, row.names=F, quote=F, sep="\t")

print("CpG Site Test")
if(is.null(dim(dmr.cols))){
	betaResults = betaRegression(formula = ~var1,
								link = "probit",
								object = predictedMeth,
								type = "BR", mc.cores=threads)
}else{
	stop("Add code for multi-variate BiSeq analysis, step #2")
}

site.fdr = p.adjust(betaResults$p.val, 'fdr')
betaResults = data.frame(betaResults[,1:3],site.fdr,betaResults[,4:ncol(betaResults)])
stat.siteID = paste(betaResults$chr,betaResults$pos,sep=":")

write.table(betaResults, dms.file, quote=F, sep="\t", row.names=F)

print("Identify Differentially Methylated Regions")
DMRs = findDMRs(betaResults,
				max.dist = max.cluster.dist,
				diff.dir = TRUE)
DMRs = data.frame(DMRs)
dmrID = paste(as.character(DMRs$seqnames),":",DMRs$start,"-",DMRs$end,sep="")


siteGR = GRanges(Rle(as.character(betaResults$chr)),
				IRanges(start=betaResults$pos, end=betaResults$pos))
regionGR = GRanges(Rle(as.character(DMRs$seqnames)),
				IRanges(start=DMRs$start, end=DMRs$end),
				Names=dmrID)
hits = findOverlaps(regionGR,siteGR)
overlaps = data.frame(pintersect(regionGR[queryHits(hits)], siteGR[subjectHits(hits)]))

sites.per.region = tapply(overlaps$Names, overlaps$Names, length)
sites.per.region = sites.per.region[match(dmrID, names(sites.per.region))]

overlap.siteID = paste(overlaps$seqnames,":",overlaps$start,sep="")
sites.in.region = tapply(overlap.siteID, overlaps$Names, paste, collapse=",")
sites.in.region = sites.in.region[match(dmrID, names(sites.in.region))]

site.fdr = betaResults$site.fdr[match(overlap.siteID,stat.siteID)]
median.fdr = tapply(site.fdr, overlaps$Names, median, na.rm=T)
median.fdr = median.fdr[match(dmrID, names(median.fdr))]

DMRs=data.frame(DMRs[1:6],median.fdr,DMRs[7:ncol(DMRs)],sites.per.region, sites.in.region)
write.table(DMRs, dmr.file, row.names=F, quote=F, sep="\t")

save.image(file=RData.file)