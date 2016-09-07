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

COHCAP.sample.table = data.frame(userID,dmr.cols)
COHCAP.sample.file = paste(paste(dmr.vars,collapse="_"),"_COHCAP_samples.txt",sep="")
write.table(COHCAP.sample.table, COHCAP.sample.file,
			col.names=F, row.names=F, sep="\t", quote=F)
			
print("Promoter DMR")
promoter.project = paste("Promoter_",project.name,sep="")
beta.table = COHCAP.annotate(beta.file, promoter.project, project.folder,
								platform="custom", annotation.file = gene.mapping,
								output.format = output.format)
								
filtered.sites = COHCAP.site(COHCAP.sample.file, beta.table, promoter.project, project.folder, ref=ref,
								methyl.cutoff=methyl.threshold, unmethyl.cutoff=unmethyl.threshold,
								delta.beta = site.delta.beta,
								pvalue.cutoff = site.pvalue, fdr.cutoff=site.fdr,
								create.wig = wig.output, paired=eval(parse(text=pair.var)),
								num.groups=num.groups, output.format = output.format)

promoter.list = COHCAP.avg.by.island(COHCAP.sample.file, filtered.sites, beta.table, promoter.project,
									project.folder, methyl.cutoff=methyl.threshold, unmethyl.cutoff = unmethyl.threshold,
									delta.beta.cutoff = island.delta.beta, pvalue.cutoff=island.pvalue, fdr.cutoff=island.fdr,
									num.groups=num.groups, num.sites=sites.per.island, plot.box=TRUE,
									paired=eval(parse(text=pair.var)), ref=ref,	output.format = output.format)						

								
print("Island DMR")
island.project = paste("CpG_Island_",project.name,sep="")
beta.table = COHCAP.annotate(beta.file, island.project, project.folder,
								platform="custom", annotation.file = island.mapping,
								output.format = output.format)
	
filtered.sites = COHCAP.site(COHCAP.sample.file, beta.table, island.project, project.folder, ref=ref,
								methyl.cutoff=methyl.threshold, unmethyl.cutoff=unmethyl.threshold,
								delta.beta = site.delta.beta,
								pvalue.cutoff = site.pvalue, fdr.cutoff=site.fdr,
								create.wig = wig.output, paired=eval(parse(text=pair.var)),
								num.groups=num.groups, output.format = output.format)

island.list = COHCAP.avg.by.island(COHCAP.sample.file, filtered.sites, beta.table, island.project,
									project.folder, methyl.cutoff=methyl.threshold, unmethyl.cutoff = unmethyl.threshold,
									delta.beta.cutoff = island.delta.beta, pvalue.cutoff=island.pvalue, fdr.cutoff=island.fdr,
									num.groups=num.groups, num.sites=sites.per.island, plot.box=TRUE,
									paired=eval(parse(text=pair.var)), ref=ref,	output.format = output.format)						
unlink(COHCAP.sample.file)
