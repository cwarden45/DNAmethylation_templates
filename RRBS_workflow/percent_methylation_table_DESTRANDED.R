min.percent.observed = 0.75

param.table = read.table("parameters.txt", header=T, sep="\t")
min.cov.total=as.numeric(as.character(param.table$Value[param.table$Parameter == "Min_Coverage"]))
min.cov.minority=as.numeric(as.character(param.table$Value[param.table$Parameter == "Min_Coverage_Pair"]))
alignment.folder=as.character(param.table$Value[param.table$Parameter == "Alignment_Folder"])
result.folder=as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
alignment.stat.file=as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
quant.type=as.character(param.table$Value[param.table$Parameter == "Quantification_Method"])
output.prefix=as.character(param.table$Value[param.table$Parameter == "methyl_percent_prefix"])
sample.description.file=as.character(param.table$Value[param.table$Parameter == "sample_description_file"])
genome=as.character(param.table$Value[param.table$Parameter == "genome"])
refFA=as.character(param.table$Value[param.table$Parameter == "FA_Ref"])

CpG.bed = paste(genome,"_CpG.bed",sep="")

if(!file.exists(CpG.bed)){
	print("Creating CpG Table")
	command = paste("perl find_CpG_sites.pl ",refFA," ",CpG.bed,sep="")
	system(command)
}#end if(!file.exists(CpG.bed))

COHCAP.methyl.file = paste(output.prefix,"_",quant.type,"_",min.cov.total,"x_for_COHCAP.txt",sep="")
stat.file = paste("summary_stats_",quant.type,".txt",sep="")
shared.sites.line.plot = paste("common_sites_per_sample_",quant.type,"_",min.cov.total,"x.png",sep="")

sep.cov.folder = paste(result.folder,"/Bismark_Coverage_Files",sep="")
dir.create(sep.cov.folder)

sample.description.table = read.table(sample.description.file, header=T, sep="\t")
alignment.stat.table = read.table(alignment.stat.file, header=T, sep="\t")
cov.files = as.character(alignment.stat.table$Cov.File[match(sample.description.table$sampleID, alignment.stat.table$Sample)])

CpG.table = read.table(CpG.bed, header=F, sep="\t")

common.Fsites=paste(CpG.table$V1,":",CpG.table$V2,sep="")
common.Rsites=paste(CpG.table$V1,":",CpG.table$V3,sep="")
rm(CpG.table)

site.count.1x=c()
site.count.4x=c()
site.count.10x=c()
site.count.20x=c()
site.count.40x=c()

if(quant.type == "Bismark"){
	for (i in 1:nrow(sample.description.table)){
		print(as.character(sample.description.table$userID[i]))
		print(cov.files[i])

		if(!file.exists(cov.files[i])){
			gzFile = paste(cov.files[i],".gz",sep="")
			if(!file.exists(gzFile)){
				stop(paste("Cannot find cov file (or .gz cov file) for ",cov.files[i],sep=""))
			}#end if(!file.exists(gzFile))
			
			command = paste("gunzip -c ",gzFile," > ",cov.files[i],sep="")
			system(command)
		}#end if(!file.exists(cov.files[i]))
		
		cov.table = read.table(cov.files[i], head=F, sep="\t")
		temp.siteID = paste(as.character(cov.table$V1),as.character(cov.table$V2),sep=":")
		temp.site.cov = cov.table$V5 + cov.table$V6
		Fcov = temp.site.cov[match(common.Fsites, temp.siteID)]
		Rcov = temp.site.cov[match(common.Rsites, temp.siteID)]
		temp.min.cov = apply(data.frame(Fcov, Rcov), 1, min)
		temp.total.cov = Fcov + Rcov
		temp.total.cov[temp.total.cov < min.cov.total]=NA
		temp.total.cov[(temp.min.cov < min.cov.minority)&!is.na(min.cov.minority)]=NA
		
		rm(temp.site.cov)
		
		Fmeth=cov.table$V5[match(common.Fsites, temp.siteID)]
		Rmeth=cov.table$V5[match(common.Rsites, temp.siteID)]
		temp.total.meth = Fmeth + Rmeth
		temp.percent = temp.total.meth / temp.total.cov
		rm(temp.total.meth)
		rm(Fcov)
		rm(Rcov)
		rm(Fmeth)
		rm(Rmeth)
		rm(temp.min.cov)
		rm(cov.table)
		
		#test.mat = data.frame(temp.total.meth, temp.total.cov, temp.percent)
		
		site.count.1x[i]=length(temp.total.cov[!is.na(temp.total.cov)])
		print(site.count.1x[i])
		site.count.4x[i]=length(temp.total.cov[!is.na(temp.total.cov)&(temp.total.cov >= 4)])
		site.count.10x[i]=length(temp.total.cov[!is.na(temp.total.cov)&(temp.total.cov >= 10)])
		site.count.20x[i]=length(temp.total.cov[!is.na(temp.total.cov)&(temp.total.cov >= 20)])
		site.count.40x[i]=length(temp.total.cov[!is.na(temp.total.cov)&(temp.total.cov >= 40)])

		rm(temp.total.cov)
		
		if(i == 1){
			percent.table = data.frame(percent=temp.percent)
			colnames(percent.table)=as.character(sample.description.table$userID[i])
			print(head(percent.table))
		}else{		
			prev.ids = colnames(percent.table)
			print(prev.ids)
			
			percent.table = data.frame(percent.table, temp.percent)
			col.ids =c(prev.ids,as.character(sample.description.table$userID[i]))
			colnames(percent.table) = col.ids
			print(head(percent.table))
		}
		
		rm(temp.percent)
		rm(temp.siteID)
		
		#comment to skip moving extracted file to result folder (if it is there already)
		command = paste("mv ",cov.files[i]," ",sep.cov.folder,sep="")
		system(command)
	}#end for (i in 1:nrow(sample.description.table))
}else if (quant.type == "methylKit"){
	stop("Need to add code to extract destranded methylation from methylKit files")
}else{
	stop(paste(quant.type," not supported type of 'Quantification_Method' (should be 'Bismark' or 'methylKit')",sep=""))
}

count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
}#end count.covered

site.covered.by.sample = apply(percent.table, 1, count.covered)
#site.covered.by.sample[site.covered.by.sample == 0]=NA
site.covered.by.sample.table = table(site.covered.by.sample)
num.samples = as.numeric(names(site.covered.by.sample.table))
site.covered.by.sample.table = site.covered.by.sample.table[order(num.samples, decreasing=T)]
num.samples = as.numeric(names(site.covered.by.sample.table))
cum.cov = c()
total.sites=0
for(i in 1:length(site.covered.by.sample.table)){
	total.sites = total.sites + site.covered.by.sample.table[i]
	cum.cov[i] = total.sites
}
cum.cov = cum.cov/1000000

png(shared.sites.line.plot)
plot(num.samples,cum.cov,
		xlab="Number of Samples", ylab="Cumulative Number of Sites (in Millions)",
		type="l", col="blue", ylim=c(0,max(cum.cov)),
		main="Commonly Covered CpGs Per Sample")
dev.off()

stat.out.table = data.frame(userID=sample.description.table$userID,
							alignment.stat.table[match(sample.description.table$sampleID,alignment.stat.table$Sample),1:(ncol(alignment.stat.table)-1)],
							site.count.1x, site.count.4x, site.count.10x, site.count.20x, site.count.40x)
write.table(stat.out.table, stat.file, quote=F, sep="\t", row.names=F)

print(dim(percent.table))
common.Fsites = common.Fsites[site.covered.by.sample > min.percent.observed * ncol(percent.table)]
percent.table = percent.table[site.covered.by.sample > min.percent.observed * ncol(percent.table),]
percent.table=round(percent.table, digits=2)
print(dim(percent.table))
percent.table = data.frame(SiteID=common.Fsites, percent.table)
write.table(percent.table, COHCAP.methyl.file, quote=F, sep="\t", row.names=F)
