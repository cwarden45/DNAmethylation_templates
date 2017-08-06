region.name = "Promoter"
region.bed = "/path/to/TxDb_hg19_promoter_1500bp.bed"
output.prefix = "TxDb_hg19_PROMOTER_1500bp_average_methylation_Bismark_10x"

#region.name = "Promoter CpG Island"
#region.bed = "/path/to/TxDb_hg19_UCSC_CpG_Islands_in_promoter_1500bp.bed"
#output.prefix = "TxDb_hg19_CpG_ISLAND_IN_PROMOTER_1500bp_average_methylation_Bismark_10x"

site.methyl.file = "../Result/percent_methylation_Bismark_10x.txt"
min.sites = 4
region.methyl.output.file = paste(output.prefix,".txt",sep="")
cov.by.sample.plot = paste(output.prefix,"_covered_samples.png",sep="")

library("GenomicRanges")

region.table = read.table(region.bed, head=F, sep="\t")
print(dim(region.table))
if(length(grep("_",region.table$V1)) > 0){
	region.table = region.table[-grep("_",region.table$V1),]
	region.table$V1 = as.factor(as.character(region.table$V1))
	print(dim(region.table))
}
refGR = GRanges(Rle(region.table$V1),
				IRanges(start=region.table$V2, end=region.table$V3),
				Names=region.table$V4,
				Rle(strand(region.table$V6)))
refGR=unique(refGR)
GR.table = data.frame(refGR)
GR.ID = paste(GR.table$seqnames,":",GR.table$start,"-",GR.table$end,sep="")

total.regions = levels(as.factor(GR.ID))

average.by.region = function(methyl.arr, region.map, min.sites){
	region.count = tapply(methyl.arr, region.map,length)
	region.mean = tapply(methyl.arr, region.map,mean, na.rm=T)
	region.mean[region.count < min.sites]=NA
	return(region.mean)
}

methyl.table = read.table(site.methyl.file, head=T, sep="\t")
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

overlap.siteID = paste(overlaps$seqnames,overlaps$start,sep=":")
regionID = as.character(overlaps$Names)

kept.sites = methyl.table$SiteID[match(overlap.siteID,methyl.table$SiteID,nomatch=0)]

print(dim(methyl.table))
methyl.table = methyl.table[match(kept.sites, as.character(methyl.table$SiteID),nomatch=0),]
regionID = regionID[match(kept.sites, overlap.siteID,nomatch=0)]
print(dim(methyl.table))

methyl.mat = methyl.table[,2:ncol(methyl.table)]
rm(methyl.table)

region.prom.methyl = apply(methyl.mat, 2, average.by.region, region.map = regionID, min.sites=min.sites)
region.prom.methyl = round(region.prom.methyl, digits=2)

region.prom.table = data.frame(Region = rownames(region.prom.methyl), region.prom.methyl)
write.table(region.prom.table, region.methyl.output.file, row.names=F, sep="\t")

#create plot of covered samples
count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
}#end count.covered

site.covered.by.sample = apply(region.prom.methyl, 1, count.covered)
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
cum.cov = cum.cov

png(cov.by.sample.plot)
plot(num.samples,cum.cov,
		xlab="Number of Samples", ylab=paste("Number of ",region.name,"s Covered",sep=""),
		type="l", col="blue", ylim=c(0,max(length(total.regions))),
		main=paste("Commonly Covered ",region.name,"s (total = ",length(total.regions),")",sep=""))
abline(h=length(total.regions),col="red")
dev.off()