param.table = read.table("parameters.txt", header=T, sep="\t")
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

project.folder = paste(overall.result.folder,"/",project.folder,sep="")

library(COHCAP)

percent.methyl.file = paste(output.prefix,"_",quant.type,"_",min.cov,"x_for_COHCAP.txt",sep="")
dir.create(project.folder)

if(!file.exists(COHCAP.mapping)){
	#create COHCAP custom mapping table
	library("GenomicRanges")
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
					
	methyl.table = read.table(percent.methyl.file, head=T, sep="\t")

	extract.value = function(char, char.index){
		char.info = unlist(strsplit(char,split=":"))
		return(char.info[char.index])
	}#end def extract.id

	site.chr=as.character(unlist(sapply(as.character(methyl.table$SiteID), extract.value, char.index=1)))
	site.pos=as.numeric(as.character(unlist(sapply(as.character(methyl.table$SiteID), extract.value, char.index=2))))

	#code based upon example from https://support.bioconductor.org/p/72656/
	testGR = GRanges(Rle(site.chr),
					IRanges(start=site.pos, end=site.pos))
	hits = findOverlaps(refGR, testGR)
	overlaps = data.frame(pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)]))
	print(dim(overlaps))

	extract.gene.region = function(char){
		char.info = unlist(strsplit(char,split="_"))
		gene = char.info[1]
		region = gsub(paste(gene,"_",sep=""),"",char)
		list.obj = list(gene,region)
		return(list.obj)
	}#end def extract.gene.region
	
	overlap.siteID = paste(overlaps$seqnames,overlaps$start,sep=":")
	overlap.table = t(sapply(as.character(overlaps$Names),extract.gene.region))
	colnames(overlap.table)=c("Gene","Island")
	
	merge.annotations = function(anno.arr){
		anno.arr = unique(as.character(anno.arr))
		if(length(anno.arr)==1){
			return(anno.arr)
		}else{
			return(paste(anno.arr,collapse=";"))
		}
		
	}#end def merge.annotations
	
	Gene = tapply(overlap.table[,"Gene"], overlap.siteID, merge.annotations)
	Gene = Gene[match(as.character(methyl.table$SiteID),names(Gene))]
	Island = tapply(overlap.table[,"Island"], overlap.siteID, merge.annotations)
	Island = Island[match(as.character(methyl.table$SiteID),names(Island))]
	
	#use similar position format as 450k/EPIC arrays
	#site.chr = gsub("chrM","MT",site.chr)
	site.chr = gsub("chr","",site.chr)
	tmp.COHCAP.region.map = data.frame(SiteID=as.character(methyl.table$SiteID),
										Chr=site.chr,Loc=site.pos,
										Gene,	Island)
	write.table(tmp.COHCAP.region.map,COHCAP.mapping,quote=F, sep="\t", row.names=F)

	#use object.size() to help determine variables that should be removed
	rm(promoter.table)
	rm(methyl.table)
	rm(tmp.COHCAP.region.map)
	rm(overlaps)
	rm(overlap.table)
	rm(Gene)
	rm(Island)
	rm(site.pos)
	rm(site.chr)
}#end if(!file.exists(COHCAP.mapping))

#other steps similar to EPIC pipeline
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

if(pair.var != "continuous"){
	pair.var = eval(parse(text=pair.var))
}	

print("Custom Region Annotations")
promoter.project = paste("Promoter_",project.name,sep="")
beta.table = COHCAP.annotate(percent.methyl.file, promoter.project, project.folder,
								platform="custom", annotation.file = COHCAP.mapping,
								output.format = output.format)

print("Custom Differentially Methylated Sites")
filtered.sites = COHCAP.site(COHCAP.sample.file, beta.table, promoter.project, project.folder, ref=ref,
								methyl.cutoff=methyl.threshold, unmethyl.cutoff=unmethyl.threshold,
								delta.beta = site.delta.beta, ttest.sub="ANOVA",
								pvalue.cutoff = site.pvalue, fdr.cutoff=site.fdr,
								create.wig = wig.output, paired=pair.var,
								num.groups=num.groups, output.format = output.format)

print("Custom Differentially Methylated Regions")
promoter.list = COHCAP.avg.by.island(COHCAP.sample.file, filtered.sites, beta.table, promoter.project,
									project.folder, methyl.cutoff=methyl.threshold, unmethyl.cutoff = unmethyl.threshold,
									delta.beta.cutoff = island.delta.beta, pvalue.cutoff=island.pvalue, fdr.cutoff=island.fdr,
									num.groups=num.groups, num.sites=sites.per.island, plot.box=TRUE, plot.heatmap =TRUE,
				     					max.cluster.dist=max.cluster.dist, ttest.sub="ANOVA",
									paired=pair.var, ref=ref, output.format = output.format)
if(expression.file != "NULL"){
	print("Custom Integration with Gene Expression Data")
	COHCAP.integrate.avg.by.island(island.list=promoter.list, project.name=promoter.project,
							project.folder=project.folder, expr.file=expression.file,
							sample.file=COHCAP.sample.file, cor.pvalue.cutoff=0.05,
							cor.fdr.cutoff=0.05, cor.cutoff=-0.2, plot.scatter=T,
							output.format = output.format, ref = ref)												
	
}#end if(expression.file != "NULL")

file.remove(COHCAP.sample.file)
