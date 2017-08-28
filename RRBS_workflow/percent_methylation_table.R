min.percent.observed = 0.75

param.table = read.table("parameters.txt", header=T, sep="\t")
min.cov=as.numeric(as.character(param.table$Value[param.table$Parameter == "Min_Coverage"]))
alignment.folder=as.character(param.table$Value[param.table$Parameter == "Alignment_Folder"])
result.folder=as.character(param.table$Value[param.table$Parameter == "Result_Folder"])
alignment.stat.file=as.character(param.table$Value[param.table$Parameter == "aligned_stats_file"])
quant.type=as.character(param.table$Value[param.table$Parameter == "Quantification_Method"])
output.prefix=as.character(param.table$Value[param.table$Parameter == "methyl_percent_prefix"])
sample.description.file=as.character(param.table$Value[param.table$Parameter == "sample_description_file"])

percent.methyl.file = paste(result.folder,"/",output.prefix,"_",quant.type,"_",min.cov,"x.txt",sep="")
COHCAP.methyl.file = paste(output.prefix,"_",quant.type,"_",min.cov,"x_for_COHCAP.txt",sep="")
stat.file = paste("summary_stats_",quant.type,".txt",sep="")
shared.sites.line.plot = paste("common_sites_per_sample_",quant.type,"_",min.cov,"x.png",sep="")

sep.cov.folder = paste(result.folder,"/Bismark_Coverage_Files",sep="")
dir.create(sep.cov.folder)

sample.description.table = read.table(sample.description.file, head=T, sep="\t")
alignment.stat.table = read.table(alignment.stat.file, head=T, sep="\t")
cov.files = as.character(alignment.stat.table$Cov.File[match(sample.description.table$sampleID, alignment.stat.table$Sample)])

common.sites=c()

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
		temp.total.cov = cov.table$V5 + cov.table$V6
		
		##use pipline="bismarkCoverage" in methylKit to avoid needing to create a re-formated coverage file (although you can uncomment code if you wish to use it for other programs)
		#methylKit.file = paste(alignment.folder,"/",sample.description.table$sampleID[i],
		#						"/",sample.description.table$sampleID[i],"_Bismark_cov_reformat.txt",sep="")
		#chrBase=paste(as.character(cov.table$V1),as.character(cov.table$V2),sep=":")
		#strand=rep("*",length(chrBase))
		#methylKit.table = data.frame(chrBase, chr=as.character(cov.table$V1),
		#							base=as.character(cov.table$V2),
		#							strand=strand,	coverage=temp.total.cov,
		#							freqC=cov.table$V4,	freqT=100-cov.table$V4)
		#write.table(methylKit.table, methylKit.file, quote=F, sep="\t", row.names=F)
		
		site.count.1x[i]=length(temp.total.cov)
		site.count.4x[i]=length(temp.total.cov[temp.total.cov >= 4])
		site.count.10x[i]=length(temp.total.cov[temp.total.cov >= 10])
		site.count.20x[i]=length(temp.total.cov[temp.total.cov >= 20])
		site.count.40x[i]=length(temp.total.cov[temp.total.cov >= 40])
		
		cov.table = cov.table[temp.total.cov >= min.cov,]
		temp.siteID = paste(as.character(cov.table$V1),as.character(cov.table$V2),sep=":")
		
		if(i == 1){
			common.sites = temp.siteID
			percent.table = data.frame(cov.table$V4)
			colnames(percent.table)=as.character(sample.description.table$userID[i])
			print(head(percent.table))
		}else{
			prev.common.sites = common.sites
			common.sites = union(common.sites, temp.siteID)
			prev.ids = colnames(percent.table)
			print(prev.ids)
			percent.table = percent.table[match(common.sites,prev.common.sites),]
			new.percent =cov.table$V4[match(common.sites, temp.siteID)]
			
			percent.table = data.frame(percent.table, new.percent)
			col.ids =c(prev.ids,as.character(sample.description.table$userID[i]))
			colnames(percent.table) = col.ids
			print(head(percent.table))
		}
		
		#move extracted file to result folder
		command = paste("mv ",cov.files[i]," ",sep.cov.folder,sep="")
		system(command)
	}#end for (i in 1:nrow(sample.description.table))
}else if (quant.type == "methylKit"){
	stop("Need to add code to extract code from methylKit files")
}else{
	stop(paste(quant.type," not supported type of 'Quantification_Method' (should be 'Bismark' or 'methylKit')",sep=""))
}

count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
}#end count.covered

site.covered.by.sample = apply(percent.table, 1, count.covered)
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

percent.table = data.frame(SiteID=common.sites, percent.table)
write.table(percent.table, percent.methyl.file, quote=F, sep="\t", row.names=F)

stat.out.table = data.frame(userID=sample.description.table$userID,
							alignment.stat.table[match(sample.description.table$sampleID,alignment.stat.table$Sample),1:(ncol(alignment.stat.table)-1)],
							site.count.1x, site.count.4x, site.count.10x, site.count.20x, site.count.40x)
write.table(stat.out.table, stat.file, quote=F, sep="\t", row.names=F)

#extra formatting for COHCAP
percent.table = percent.table[,2:ncol(percent.table)]
print(dim(percent.table))
common.sites = common.sites[site.covered.by.sample > min.percent.observed * ncol(percent.table)]
percent.table = percent.table[site.covered.by.sample > min.percent.observed * ncol(percent.table),]
print(dim(percent.table))
percent.table = data.frame(SiteID=common.sites, percent.table)
write.table(percent.table, COHCAP.methyl.file, quote=F, sep="\t", row.names=F)
