param.table = read.table("parameters.txt", header=T, sep="\t")
sample.description.file = as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
beta.prefix  = as.character(param.table$Value[param.table$Parameter == "beta_prefix"])
beta.normalization  = as.character(param.table$Value[param.table$Parameter == "beta_normalization"])
dmr.vars  = as.character(param.table$Value[param.table$Parameter == "dmr_groups"])
output.format  = as.character(param.table$Value[param.table$Parameter == "COHCAP_output_format"])
project.folder = as.character(param.table$Value[param.table$Parameter == "project_folder"])
project.name = as.character(param.table$Value[param.table$Parameter == "comp_name"])
expression.file = as.character(param.table$Value[param.table$Parameter == "expression_file"])
island.mapping = as.character(param.table$Value[param.table$Parameter == "island_mapping"])
gene.mapping = as.character(param.table$Value[param.table$Parameter == "gene_mapping"])
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

library(COHCAP)

beta.file = paste(beta.prefix,"_",beta.normalization,".txt",sep="")
dir.create(project.folder)

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

print("Promoter to 1st Exon Annotations")
promoter.project = paste("Promoter_",project.name,sep="")
beta.table = COHCAP.annotate(beta.file, promoter.project, project.folder,
								platform="custom", annotation.file = gene.mapping,
								output.format = output.format)

print("Promoter to 1st Exon Differentially Methylated Sites")
filtered.sites = COHCAP.site(COHCAP.sample.file, beta.table, promoter.project, project.folder, ref=ref,
								methyl.cutoff=methyl.threshold, unmethyl.cutoff=unmethyl.threshold,
								delta.beta = site.delta.beta,
								pvalue.cutoff = site.pvalue, fdr.cutoff=site.fdr,
								create.wig = wig.output, paired=pair.var,
								num.groups=num.groups, output.format = output.format)

print("Promoter to 1st Exon Differentially Methylated Regions")
promoter.list = COHCAP.avg.by.island(COHCAP.sample.file, filtered.sites, beta.table, promoter.project,
									project.folder, methyl.cutoff=methyl.threshold, unmethyl.cutoff = unmethyl.threshold,
									delta.beta.cutoff = island.delta.beta, pvalue.cutoff=island.pvalue, fdr.cutoff=island.fdr,
									num.groups=num.groups, num.sites=sites.per.island, plot.box=TRUE, plot.heatmap =TRUE,
				     					max.cluster.dist=max.cluster.dist,
									paired=pair.var, ref=ref, output.format = output.format)
if(expression.file != "NULL"){
	print("Promoter to 1st Exon Integration with Gene Expression Data")
	COHCAP.integrate.avg.by.island(island.list=promoter.list, project.name=promoter.project,
							project.folder=project.folder, expr.file=expression.file,
							sample.file=COHCAP.sample.file, cor.pvalue.cutoff=0.05,
							cor.fdr.cutoff=0.05, cor.cutoff=-0.2, plot.scatter=T,
							output.format = output.format, ref = ref)												
	
}#end if(expression.file != "NULL")

#promoter results likely sufficient
stop()	

print("UCSC Island Annotations")
island.project = paste("CpG_Island_",project.name,sep="")
beta.table = COHCAP.annotate(beta.file, island.project, project.folder,
								platform="custom", annotation.file = island.mapping,
								output.format = output.format)
	
filtered.sites = COHCAP.site(COHCAP.sample.file, beta.table, island.project, project.folder, ref=ref,
								methyl.cutoff=methyl.threshold, unmethyl.cutoff=unmethyl.threshold,
								delta.beta = site.delta.beta,
								pvalue.cutoff = site.pvalue, fdr.cutoff=site.fdr,
								create.wig = wig.output, paired=pair.var,
								num.groups=num.groups, output.format = output.format)

island.list = COHCAP.avg.by.island(COHCAP.sample.file, filtered.sites, beta.table, island.project,
									project.folder, methyl.cutoff=methyl.threshold, unmethyl.cutoff = unmethyl.threshold,
									delta.beta.cutoff = island.delta.beta, pvalue.cutoff=island.pvalue, fdr.cutoff=island.fdr,
									num.groups=num.groups, num.sites=sites.per.island, plot.box=TRUE, plot.heatmap =TRUE,
				   					max.cluster.dist=max.cluster.dist,
									paired=pair.var, ref=ref, output.format = output.format)						
unlink(COHCAP.sample.file)
