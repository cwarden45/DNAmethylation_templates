comp.name = "avg_GENE_[comp name]"
percent.file = "TxDb_hg19_GENE_promoter_1500bp_average_methylation_Bismark_10x.txt"
extract.gene.flag=FALSE

#comp.name = "avg_CPG_ISLAND_[comp name]"
#extract.gene.flag=TRUE
#percent.file = "TxDb_hg19_CpG_ISLAND_IN_PROMOTER_1500bp_average_methylation_Bismark_10x.txt"

sample.description.file = "sample_description.txt"
trt.group = "Trt"
plot.groups = c("Group")
dmg.groups = c("Group")

#p-value method can be 'ANOVA', 'lm' (linear regression), or 'limma'
pvalue.method = "limma"
dp.cutoff = 20
pvalue.cutoff = 0.05
fdr.cutoff = 0.25

### ideally, edit above this point ###
plot.min=0
plot.max=100
cluster.distance = "Pearson_Dissimilarity"
min.fraction.covered = 0.75
fdr.method = "q-value"
interaction.flag = "no"
lower.cont.quantile=0
upper.cont.quantile=1

avgGroupMethylation = function (geneMethyl, groups) {
	avg.methyl = tapply(geneMethyl, groups, mean, na.rm=T)
	return(avg.methyl)
}#end def avgGroupMethylation

count.defined.values = function(arr)
{
	return(length(arr[!is.na(arr)]))
}#end def count.values


gene.lm = function(arr, var1, var2=c(), var3=c())
{	
	if (length(var2) == 0){
		fit = lm(as.numeric(arr) ~ var1)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else if (length(var3) == 0){
		fit = lm(as.numeric(arr) ~ var1 + var2)
		result = summary(fit)
		pvalue = result$coefficients[2,4]
	} else {
		fit = lm(as.numeric(arr) ~ var2*var3 + var2 + var3)
		result = summary(fit)
		pvalue = result$coefficients[4,4]
	}
	return(pvalue)
}#end def gene.lm

gene.aov = function(arr, var1, var2=c(), var3=c())
{	
	#probably need to update for other variables (or require coverage in all samples)
	measured.values = var1[!is.na(arr)]
	grp.no.na = as.factor(as.character(var1[!is.na(arr)]))
	if((length(levels(grp.no.na)) < 2)|(length(measured.values) < 4)){
		return(1)
	}else{
		if (length(var2) == 0){
			fit = aov(as.numeric(arr) ~ var1)
			result = summary(fit)
			aov.pvalue = result[[1]][['Pr(>F)']][1]
		} else if (length(var3) == 0){
			fit = aov(as.numeric(arr) ~ var1 + var2)
			result = summary(fit)
			aov.pvalue = result[[1]][['Pr(>F)']][1]
		} else {
			fit = aov(as.numeric(arr) ~ var2*var3 + var2 + var3)
			result = summary(fit)
			aov.pvalue = result[[1]][['Pr(>F)']][3]
		}
		return(aov.pvalue)
	}
}#end def gene.aov

calc.gene.cor = function(arr, indep.var)
{	
	na.count = length(arr[!is.na(arr)])
	if((na.count >= 3) & (sd(arr) != 0)){
		gene.cor.coef = cor(arr,indep.var)
	} else {
		gene.cor.coef = NA
	}
	return(gene.cor.coef)
}#end def calc.gene.cor

cor.dist = function(mat){
	cor.mat = cor(as.matrix(t(mat)))
	dis.mat = 1 - cor.mat
	return(as.dist(dis.mat))	
}#end def cor.dist

extract.gene.from.region = function(char){
	char.info = unlist(strsplit(char,split="_"))
	return(char.info[1])
}#end def extract.gene.from.region

fixed.color.palatte = c("green","orange","purple","cyan","pink","maroon","yellow","grey","red","blue","black",rainbow(100))
continuous.color.breaks = 10

sample.description.table = read.table(sample.description.file, sep="\t", header=T)
userID = as.character(sample.description.table$userID)

deg.group.table = sample.description.table[,dmg.groups]
if (length(dmg.groups) == 1){
	deg.meta = sample.description.table[!is.na(deg.group.table),]
} else {
	deg.grp.na.counts = apply(deg.group.table, 1, count.na.values)
	deg.meta = sample.description.table[deg.grp.na.counts == 0,]
}

percent.table = read.table(percent.file, sep="\t", header=T)
regionID = percent.table[,1]
percent.mat = percent.table[,match(userID, names(percent.table))]

genes.percent.detected = 100*apply(percent.mat, 1, count.defined.values) / ncol(percent.mat)
print(dim(percent.mat))
percent.mat = percent.mat[genes.percent.detected >= min.fraction.covered,]
regionID = regionID[genes.percent.detected >= min.fraction.covered]
print(dim(percent.mat))
rownames(percent.mat) = regionID

if(extract.gene.flag){
	gene=unlist(sapply(as.character(regionID),extract.gene.from.region))
}else{
	gene = regionID
}
if(length(plot.groups) == 1){
	print("Averaging Percent Methylation for One Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups]
} else if ((length(plot.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Percent Methylation for First Variable (for plot.groups)")
	grp = sample.description.table[,plot.groups[1]]
} else if (length(plot.groups) == 2){
	print("Averaging Percent Methylation for Interaction Variable (for plot.groups)")
	stop("Double-check code for interaction")
	grp = paste(sample.description.table[,plot.groups[1]],sample.description.table[,plot.groups[2]],sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

groupIDs = as.character(levels(as.factor(grp)))
average.percent = data.frame(t(apply(percent.mat, 1, avgGroupMethylation, groups = grp)))
if(length(groupIDs) == 1){
	average.percent = t(average.percent)
} else {
	average.percent = average.percent
}
colnames(average.percent) = paste("avg.percent", sub("-",".",groupIDs), sep=".")

#remove undefined group IDs (so, you can visualize samples that aren't in your comparison)
if(length(dmg.groups) == 1){
	var1 = sample.description.table[,dmg.groups]
	dmg.percent = percent.mat[,!is.na(var1)]
	var1 = var1[!is.na(var1)]
	if (trt.group != "continuous"){
		var1 = as.factor(as.character(var1[!is.na(var1)]))
	}
} else if (length(dmg.groups) == 2){
		var1 = sample.description.table[,dmg.groups[1]]
		var2 = sample.description.table[,dmg.groups[2]]
		dmg.samples = !is.na(var1)&!is.na(var2)
		dmg.percent = percent[,dmg.samples]
		var1 = var1[dmg.samples]
		if (trt.group != "continuous"){
			var1 = as.factor(as.character(var1[!is.na(var1)]))
		}
		var2 = var2[dmg.samples]
		if (trt.group2 != "continuous"){
			var2 = as.factor(as.character(var2[!is.na(var2)]))
		}
} else {
	stop("Code currently doesn't support more than 2 group model for DMG (with or without interaction)")
}

if(length(dmg.groups) == 1){
	print("Averaging Percent Methylation for One Variable (for dmg.groups)")
	contrast.grp = var1
} else if ((length(dmg.groups) == 2)&(interaction.flag == "no")){
	print("Averaging Percent Methylation for First Variable (for dmg.groups)")
	contrast.grp = var1
} else if (length(dmg.groups) == 2){
	print("Averaging Percent Methylation for Interaction Variable (for dmg.groups)")
	contrast.grp = paste(var1,var2,sep=":")
} else {
	stop("Code only compatible with 2 variables (with or without a 3rd interaction variable")
}

if (trt.group == "continuous"){
	contrast.grp = as.numeric(contrast.grp)
	
	gene.cor = apply(dmg.percent, 1, calc.gene.cor, indep.var=contrast.grp)

	percent.min= apply(dmg.percent, 1, quantile, na.rm=TRUE, probs=c(lower.cont.quantile))
	percent.max= apply(dmg.percent, 1, quantile, na.rm=TRUE, probs=c(upper.cont.quantile))
			
	upper.percent = percent.max
	upper.percent[!is.na(percent.cor) & (percent.cor < 0)] = percent.min[!is.na(percent.cor) & (percent.cor < 0)]
	lower.percent = percent.min
	lower.percent[!is.na(percent.cor) & (percent.cor < 0)] = percent.max[!is.na(percent.cor) & (percent.cor < 0)]
			
	rm(percent.min)
	rm(percent.max)
			
	dp = upper.percent - lower.percent
	
	fc.table = data.frame(cor=gene.cor, delta.percent=dp)
} else {
	groupIDs = as.character(levels(as.factor(contrast.grp)))
	contrast.percent = data.frame(t(apply(dmg.percent, 1, avgGroupMethylation, groups = contrast.grp)))
	colnames(contrast.percent) = paste("avg.percent", sub("-",".",groupIDs), sep=".")
}#end else

if((interaction.flag == "no") & (trt.group != "continuous")){
	print("Calculating delta-percent for primary variable")
	trt.percent = contrast.percent[,paste("avg.percent", sub("-",".",trt.group), sep=".")]
	cntl.percent = contrast.percent[,paste("avg.percent", sub("-",".",groupIDs[groupIDs != trt.group]), sep=".")]

	dp = round(trt.percent - cntl.percent, digits = 4)
	fc.table = data.frame(delta.percent=dp)
} else{
		stop("Add extra code?")
}#end else

rep.check = 1
for (i in 1:length(dmg.groups)){
	deg.group = dmg.groups[i]
	
	if((i == 1) & (trt.group != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if ((i == 2) & (trt.group2 != "continuous")){
		deg.group.values = as.factor(as.character(deg.meta[,deg.group]))
		min.reps = min(table(deg.group.values))
		if (min.reps < 2){
			rep.check=0
			print("There are not at least 2 samples per-group in order to calculate p-value.")
			print("In the future, please make sure you at least have duplicate samples.")
		}#end if (min.reps < 2)
	} else if (i > 2){
		stop("Workflow currently doesn't support use of more than 2 variables")
	}
}#end for (deg.group in dmg.groups)

if(rep.check == 1){
	#start p-value calculation
	if (pvalue.method == "limma"){
		library(limma)

		if (length(dmg.groups) == 1){
			print("limma with 1 variable")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			}
			design = model.matrix(~var1)
			fit = lmFit(dmg.percent,design)
			fit = eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			test.pvalue = pvalue.mat[,2]
		} else if ((length(dmg.groups) == 2)&(interaction.flag == "no")){
			print("limma with 2 variables")

			if (trt.group == "continuous"){
				var1 = as.numeric(var1)
			} else{
				var1 = as.factor(as.character(var1))
			}

			if (trt.group2 == "continuous"){
				var2 = as.numeric(var2)
			} else{
				var2 = as.factor(as.character(var2))
			}
			design = model.matrix(~var1 + var2)
			fit = lmFit(dmg.percent,design)
			fit = eBayes(fit)
			pvalue.mat = data.frame(fit$p.value)
			test.pvalue = pvalue.mat[,2]
		}#end else if ((length(dmg.groups) == 2)&(interaction.flag == "no"))
	} else if (pvalue.method == "lm"){
			if (length(dmg.groups) == 1){
				print("linear regression with 1 variable")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				test.pvalue = apply(dmg.percent, 1, gene.lm, var1=var1)
			} else if ((length(dmg.groups) == 2)&(interaction.flag == "no")){
				print("linear regression with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(dmg.percent, 1, gene.lm, var1=var1, var2=var2)
			} else if ((length(dmg.groups) == 2)&(interaction.flag == "model")){
				print("linear regression with 2 variables plus interaction")
				var3 = as.factor(paste(var1,var2,sep=":"))

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(dmg.percent, 1, gene.lm, var1=var3, var2=var1, var3=var2)
			}
	} else if (pvalue.method == "ANOVA"){
		dmg.percent = as.matrix(dmg.percent)
			if (length(dmg.groups) == 1){
				print("ANOVA with 1 variable")
				
				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}
				test.pvalue = apply(dmg.percent, 1, gene.aov, var1=var1)
			} else if ((length(dmg.groups) == 2)&(interaction.flag == "no")){
				print("ANOVA with 2 variables")

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(dmg.percent, 1, gene.aov, var1=var1, var2=var2)
			} else if ((length(dmg.groups) == 2)&(interaction.flag == "model")){
				print("ANOVA with 2 variables plus interaction")
				var3 = as.factor(paste(var1,var2,sep=":"))

				if (trt.group == "continuous"){
					var1 = as.numeric(var1)
				}

				if (trt.group2 == "continuous"){
					var2 = as.numeric(var2)
				}
				test.pvalue = apply(dmg.percent, 1, gene.aov, var1=var3, var2=var1, var3=var2)
			}
	} else{
		stop("pvalue_method must be \"limma\", \"lm\", or \"ANOVA\"")
	}
} else{
	test.pvalue = rep(1,times=length(regionID))
	prim.pvalue = rep(1,times=length(regionID))
	sec.pvalue = rep(1,times=length(regionID))
}#end else

if (trt.group == "continuous"){
	upID = "Increased Methylation"
	downID = "Decreased Methylation"
} else {
	upID = paste(trt.group," Up",sep="")
	downID = paste(trt.group," Down",sep="")	
}


if (interaction.flag == "no"){
	if (fdr.method == "BH"){
		fdr = p.adjust(test.pvalue, "fdr")
	} else if (fdr.method == "q-value"){
		library(qvalue)
		qobj = qvalue(p = test.pvalue)
		fdr = qobj$qvalue
		png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else if (fdr.method == "q-lfdr"){
		library(qvalue)
		qobj = qvalue(p = test.pvalue)
		fdr = qobj$lfdr
		png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
		qHist = hist(qobj)
		print(qHist)
		dev.off()
	} else {
		stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
	}
	status = rep("No Change", times=length(fdr))
	if (trt.group == "continuous"){
		status[(gene.cor >= dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(gene.cor <= -dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(pvalue = test.pvalue, qvalue = fdr)
	} else{
		status[(dp >= dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
		status[(dp <= -dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
		pvalue.table = data.frame(pvalue = test.pvalue, qvalue = fdr)
	}#end else
} else{
	trt.group = prim.trt
	if(interaction.flag == "model"){
		if (fdr.method == "BH"){
			fdr = p.adjust(test.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj = qvalue(p = test.pvalue)
			fdr = qobj$qvalue
			png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj = qvalue(p = test.pvalue)
			fdr = qobj$lfdr
			png(paste(comp.name,"_",pvalue.method,"_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		status = rep("No Change", times=length(fdr))
		if ((trt.group == "continuous")&(trt.group2 == "continuous")){
			status[(gene.cor.int >= dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(gene.cor.int <= -dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(pvalue = test.pvalue, qvalue = fdr)
		} else if ((trt.group != "continuous")&(trt.group2 != "continuous")){
			status[(overall.dp >= dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			status[(overall.dp <= -dp.cutoff) & (test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = downID
			pvalue.table = data.frame(pvalue = test.pvalue, qvalue = fdr)
		} else {
			upID = "Variable Methylation"
			status[(test.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = upID
			pvalue.table = data.frame(pvalue = test.pvalue, qvalue = fdr)
		}#end else
	} else if (interaction.flag == "filter-overlap"){
		if (fdr.method == "BH"){
			fdr = p.adjust(prim.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj = qvalue(p = prim.pvalue)
			fdr = qobj$qvalue
			png(paste(comp.name,"_",pvalue.method,"_prim_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj = qvalue(p = prim.pvalue)
			fdr = qobj$lfdr
			png(paste(comp.name,"_",pvalue.method,"_prim_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}
		pass1.status = rep("No Change", times=length(fdr))
		if (trt.group == "continuous"){
			pass1.status[(prim.dp >= dp.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(prim.dp <= -dp.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		} else{
			pass1.status[(prim.dp >= dp.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Up",sep="")
			pass1.status[(prim.dp <= -dp.cutoff) & (prim.pvalue <= pvalue.cutoff) & (fdr <= fdr.cutoff)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Primary Up-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Primary Down-Regulated: ",length(pass1.status[pass1.status == paste(trt.group," Down",sep="")]),sep=""))

		if (fdr.method == "BH"){
			sec.fdr = p.adjust(sec.pvalue, "fdr")
		} else if (fdr.method == "q-value"){
			library(qvalue)
			qobj = qvalue(p = sec.pvalue)
			sec.fdr = qobj$qvalue
			png(paste(comp.name,"_",pvalue.method,"_sec_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else if (fdr.method == "q-lfdr"){
			library(qvalue)
			qobj = qvalue(p = sec.pvalue)
			sec.fdr = qobj$lfdr
			png(paste(comp.name,"_",pvalue.method,"_sec_qvalue_plot.png",sep=""))
			qHist = hist(qobj)
			print(qHist)
			dev.off()
		} else {
			stop("fdr_method must be \"BH\", \"q-value\", or \"q-lfdr\"")
		}		

		pass2.status = rep("No Change", times=length(fdr))
		if (trt.group2 == "continuous"){
			pass2.status[(gene.cor2 >= dp.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(gene.cor2 <= -dp.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		} else{
			pass2.status[(sec.dp >= dp.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Up",sep="")
			pass2.status[(sec.dp <= -dp.cutoff2) & (sec.pvalue <= pvalue.cutoff2) & (sec.fdr <= fdr.cutoff2)] = paste(trt.group," Down",sep="")
		}#end else

		print(paste("Secondary Up-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Up",sep="")]),sep=""))
		print(paste("Secondary Down-Regulated: ",length(pass2.status[pass2.status == paste(trt.group," Down",sep="")]),sep=""))
			
		pvalue.table = data.frame(prim.pvalue = prim.pvalue, prim.FDR = fdr,
										sec.pvalue=sec.pvalue, sec.fdr=sec.fdr)
			
		status = rep("No Change", times=length(fdr))
		status[(pass1.status == paste(trt.group," Up",sep="")) & (pass2.status == "No Change")] = upID
		status[(pass1.status == paste(trt.group," Down",sep="")) & (pass2.status == "No Change")] = downID
	} else{
		stop("interaction must be \"no\", \"model\", or \"filter-overlap\"")
	}#end else
}#end else

print(paste("Up-Regulated: ",length(status[status == upID]),sep=""))
print(paste("Down-Regulated: ",length(status[status == downID]),sep=""))

if (interaction.flag == "filter-overlap"){
	pvalue.method = paste(pvalue.method,"two-step_filtered",sep="_")
}

if(rep.check == 1){
	deg.table = data.frame(regionID, gene,
							average.percent, fc.table,
							pvalue.table, status = status)
} else {
	deg.table = data.frame(regionID, gene,
							average.percent, fc.table, status = status)	
}#end else

deg.file = paste(comp.name,"_",pvalue.method,"_DMG_fc_",dp.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".txt",sep="")
deg.file = gsub(":",".",deg.file)
write.table(deg.table, file=deg.file, row.names=F, quote=F, sep="\t")

temp.percent = percent.mat
temp.percent = temp.percent[status != "No Change", ]
dmg.genes = gene[status != "No Change"]

		num.breaks = 33
		plot.range = plot.max - plot.min
		heatmap.breaks = seq(plot.min, plot.max, by=plot.range/num.breaks)

if(length(dmg.genes) > 1){
	if (cluster.distance == "Pearson_Dissimilarity"){
		print("Using Pearson Dissimilarity as Distance in Heatmap...")
		dist.fun = cor.dist
	}else{
		dist.fun=dist
	}
	
	if(length(plot.groups) > 1){
		source("heatmap.3.R")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
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
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			grp1 = as.numeric(sample.description.table[,plot.groups[1]])
			grp2 = as.character(sample.description.table[,plot.groups[2]])
		
			labelColors1 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp1)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors1[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
			
			group.levels = c(levels(as.factor(grp2)))
			color.palette = fixed.color.palatte[3:(2+length(group.levels))]
			labelColors2 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors2[grp2 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}else{
			grp1 = as.character(sample.description.table[,plot.groups[1]])
			grp2 = as.numeric(sample.description.table[,plot.groups[2]])

			group.levels = c(levels(as.factor(grp1)))
			color.palette = fixed.color.palatte[1:(length(group.levels))]
			labelColors1 = rep("black",times=length(sample.label))
			for (i in 1:length(group.levels)){
				labelColors1[grp1 == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
			
			labelColors2 = rep("black",times=length(sample.label))
			library("RColorBrewer")
			continuous.color.breaks = 10
				
			plot.var = as.numeric(grp2)
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
				
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
				
			color.range = colorRampPalette(c("purple","black","cyan"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors2[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}
		
		temp.percent = t(temp.percent)
		if(length(dmg.genes) < 25){
			col.labels = dmg.genes
		} else {
			col.labels = rep("", length(dmg.genes))
		}
		rownames(temp.percent) = sample.label
		colnames(temp.percent) = rep("", length(dmg.genes))

		column_annotation = as.matrix(dmg.genes)
		colnames(column_annotation) = c("")

		row_annotation = data.frame(label1 = labelColors1, label2 = labelColors2)
		row_annotation = as.matrix(t(row_annotation))
		rownames(row_annotation) = c(plot.groups)

		heatmap.file = paste(comp.name,"_",pvalue.method,"_DEG_fc_",dp.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.3(temp.percent,   distfun = dist.fun, hclustfun = hclust,
				breaks = heatmap.breaks,labCol=col.labels,
			  col=colorpanel(33, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
			RowSideColors=row_annotation, trace="none", margins = c(10,15),RowSideColorsSize=4, dendrogram="both")
		if((trt.group != "continuous")&(trt.group2 != "continuous")){
					legend("topright", legend=group.levels,
							col=color.palette,
							pch=15, cex=0.7)
		}else if((trt.group == "continuous")&(trt.group2 == "continuous")){
			stop("Add code for two continuous variables")
		}else if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}else{
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
			legend("topright", legend=group.levels, col=color.palette, pch=15, cex=0.7)
		}
		dev.off()
			
		if(interaction.flag != "no"){
			temp.fc.table = as.matrix(fc.table)
			if (((trt.group == "continuous") & (trt.group2 == "continuous")) | ((trt.group != "continuous") & (trt.group2 != "continuous"))){
				temp.fc.table = temp.fc.table[,-ncol(temp.fc.table)]
			}
			temp.fc.table = temp.fc.table[status != "No Change", ]
			if(length(dmg.genes) < 25){
				rownames(temp.fc.table) = dmg.genes
			} else {
				rownames(temp.fc.table) = rep("",times=length(dmg.genes))
			}
			colnames(temp.fc.table) = gsub(".:.",":",gsub("fold.change.","",colnames(temp.fc.table)))
		
			temp.fc.table[temp.fc.table < -10] = -10
			temp.fc.table[temp.fc.table > 10] = 10
		
			heatmap.file = paste("fold_change_",comp.name,"_",pvalue.method,"_DEG_fc_",dp.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
			heatmap.file = gsub(":",".",heatmap.file)
			png(file = heatmap.file)
			heatmap.2(temp.fc.table,  distfun = dist.fun, hclustfun = hclust,
				  col=colorpanel(num.breaks, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
				  trace="none", margins = c(20,5), cexCol=1.5)
			dev.off()
		}#end if(interaction.flag != "no")
		
	} else {
		labelColors = rep("black",times=length(sample.label))
		if(trt.group == "continuous"){
			library("RColorBrewer")
			continuous.color.breaks = 10
			
			plot.var = as.numeric(sample.description.table[,plot.groups])
			plot.var.min = min(plot.var, na.rm=T)
			plot.var.max = max(plot.var, na.rm=T)
			
			plot.var.range = plot.var.max - plot.var.min
			plot.var.interval = plot.var.range / continuous.color.breaks
			
			color.range = colorRampPalette(c("green","black","orange"))(n = continuous.color.breaks)
			plot.var.breaks = plot.var.min + plot.var.interval*(0:continuous.color.breaks)
			for (j in 1:continuous.color.breaks){
				#print(paste(plot.var.breaks[j],"to",plot.var.breaks[j+1]))
				labelColors[(plot.var >= plot.var.breaks[j]) &(plot.var <= plot.var.breaks[j+1])] = color.range[j]
			}#end for (j in 1:continuous.color.breaks)
		}else{
			group.levels = levels(as.factor(sample.description.table[,plot.groups]))
			color.palette = fixed.color.palatte[1:length(group.levels)]
			for (i in 1:length(group.levels)){
				labelColors[grp == as.character(group.levels[i])] = color.palette[i]
			}#end for (i in 1:length(group.levels))
		}

		temp.percent = t(temp.percent)
		if(length(dmg.genes) < 25){
			col.labels = dmg.genes
		} else {
			col.labels = rep("", length(dmg.genes))
		}
		rownames(temp.percent) = sample.label
		colnames(temp.percent) = rep("", length(dmg.genes))
			
		heatmap.file = paste(comp.name,"_",pvalue.method,"_DEG_fc_",dp.cutoff,"_fdr_",fdr.cutoff,"_pval_",pvalue.cutoff,".png",sep="")
		heatmap.file = gsub(":",".",heatmap.file)
		png(file = heatmap.file)
		heatmap.2(temp.percent, distfun = dist.fun, hclustfun = hclust,
				breaks = heatmap.breaks,labCol=col.labels,
			  col=colorpanel(num.breaks, low="blue", mid="black", high="red"), density.info="none", key=TRUE,
			 RowSideColors=labelColors, trace="none", margins = c(5,15))

		if(trt.group == "continuous"){
			legend("right",legend=c(round(plot.var.max,digits=1),rep("",length(color.range)-2),round(plot.var.min,digits=1)),
								col=rev(color.range),  pch=15, y.intersp = 0.4, cex=0.8, pt.cex=1.5)
		}else{
			legend("topright", group.levels, col=color.palette, pch=15)
		}
		dev.off()
	}#end else
}#end if(length(dmg.genes) > 1)