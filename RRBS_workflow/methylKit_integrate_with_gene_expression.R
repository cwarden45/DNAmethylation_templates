min.percent.observed = 0.75
#add extra filter to avoid getting errors in samples with too few measurements

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

library("COHCAP")

project.folder = paste(overall.result.folder,"/",project.folder,sep="")

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
beta.mat = percent.table[,6:ncol(percent.table)]/100
beta.table = data.frame(island=percent.table$regionID, gene = gene,beta.mat)
print(dim(beta.table))
count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
}#end count.covered
covered.sample.count = apply(beta.mat, 1, count.covered)
beta.table = beta.table[covered.sample.count >= min.percent.observed * ncol(beta.mat),]
print(dim(beta.table))

#reformat dmr - assumes 2-group 
dmr.table = read.table(dmr.file, head=T, sep="\t")
gene=as.character(unlist(sapply(as.character(dmr.table$regionID), extract.gene.from.region)))

avg.cols = names(dmr.table)[grep(".avg.percent.methylation",names(dmr.table))]
avg.cols = gsub(".avg.percent.methylation","",avg.cols)
ref.group = avg.cols[avg.cols != trt.group]

avg.treatment = dmr.table[,paste(trt.group,".avg.percent.methylation",sep="")] / 100
avg.reference = dmr.table[,paste(ref.group,".avg.percent.methylation",sep="")] / 100
delta.beta = dmr.table$meth.diff/100

#just fake number of sites
filtered.island.stats = data.frame(island=dmr.table$regionID, gene = gene,
								avg.treatment, avg.reference, delta.beta,
								island.pvalue=dmr.table$pvalue, island.fdr=dmr.table$qvalue,
								num.sites=rep(sites.per.island,length(gene)))
print(dim(filtered.island.stats))
filtered.island.stats = filtered.island.stats[match(beta.table$island, filtered.island.stats$island,nomatch=0),]
print(dim(filtered.island.stats))
						
island.list = list(beta.table=beta.table, filtered.island.stats=filtered.island.stats)

#create COHCAP-style sample table
sample.table = read.table(sample.description.file, sep="\t", header=T)
userID = as.character(sample.table$userID)
dmr.vars = unlist(strsplit(dmr.vars,","))
dmr.cols = sample.table[,dmr.vars]

ref = ""
if(trt.group == "continuous"){
	ref = trt.group
} else{
	if(is.null(dim(dmr.cols))){
		grp.levels = levels(as.factor(dmr.cols))
		ref = grp.levels[grp.levels != trt.group]
		if ((length(ref) != 1)&(num.groups==2)){
			stop("Reconsider design: You specified 2 groups, but provide more than 2 groups")
		}
	} else{
		grp.levels = levels(as.factor(dmr.cols[[1]]))
		ref = grp.levels[grp.levels != trt.group]
		if ((length(ref) != 1)&(num.groups==2)){
			stop("Reconsider design: You specified 2 groups, but provide more than 2 groups")
		}
	}
}#end else

if(is.null(dim(dmr.cols))){
	names(dmr.cols) = userID
	dmr.cols = dmr.cols[!is.na(dmr.cols)]
	COHCAP.sample.table = data.frame(userID=names(dmr.cols),dmr.cols)
	COHCAP.sample.file = paste(project.folder,"/",project.name,"_COHCAP_samples.txt",sep="")
	write.table(COHCAP.sample.table, COHCAP.sample.file,
				col.names=F, row.names=F, sep="\t", quote=F)
}else{
	rownames(dmr.cols) = userID
	dmr.cols = na.omit(dmr.cols)
	COHCAP.sample.table = data.frame(userID=rownames(dmr.cols),dmr.cols)
	COHCAP.sample.file = paste(project.folder,"/",project.name,"_COHCAP_samples.txt",sep="")
	write.table(COHCAP.sample.table, COHCAP.sample.file,
				col.names=F, row.names=F, sep="\t", quote=F)
}

COHCAP.integrate.avg.by.island(island.list=island.list, project.name=project.name,
							project.folder=project.folder, expr.file=expression.file,
							sample.file=COHCAP.sample.file, cor.pvalue.cutoff=0.05,
							cor.fdr.cutoff=0.05, cor.cutoff=-0.2, plot.scatter=T,
							output.format = output.format)	

file.remove(COHCAP.sample.file)