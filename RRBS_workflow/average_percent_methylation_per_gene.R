min.sites = 4
site.percent.file = "../Result/percent_methylation_Bismark_10x.txt"
gene.mapping.file = "TxDb_hg19_promoter_1500bp_COHCAP_mapping.txt"
output.prefix = "TxDb_hg19_GENE_promoter_1500bp_average_methylation_Bismark_10x"

region.percent.output.file = paste(output.prefix,".txt",sep="")
cov.by.sample.plot = paste(output.prefix,"_covered_samples.png",sep="")

average.by.gene = function(percent.arr, gene.map, min.sites){
	gene.count = tapply(percent.arr, gene.map,length)
	gene.mean = tapply(percent.arr, gene.map,mean, na.rm=T)
	gene.mean[gene.count < min.sites]=NA
	return(gene.mean)
}

gene.mapping.table = read.table(gene.mapping.file, head=T, sep="\t")
total.genes = levels(gene.mapping.table$Gene)

percent.table = read.table(site.percent.file, head=T, sep="\t")
kept.sites = percent.table$SiteID[match(gene.mapping.table$SiteID,percent.table$SiteID,nomatch=0)]

print(dim(percent.table))
percent.table = percent.table[match(kept.sites, as.character(percent.table$SiteID),nomatch=0),]
gene.mapping.table = gene.mapping.table[match(kept.sites, as.character(gene.mapping.table$SiteID),nomatch=0),]
print(dim(percent.table))

percent.mat = percent.table[,2:ncol(percent.table)]
rm(percent.table)

gene.prom.percent = apply(percent.mat, 2, average.by.gene, gene.map = gene.mapping.table$Gene, min.sites=min.sites)
gene.prom.percent = round(gene.prom.percent, digits=2)

gene.prom.table = data.frame(Gene = rownames(gene.prom.percent), gene.prom.percent)
write.table(gene.prom.table, region.percent.output.file, row.names=F, sep="\t")

#create plot of covered samples
count.covered = function(arr){
	return(length(arr[!is.na(arr)]))
}#end count.covered

site.covered.by.sample = apply(gene.prom.percent, 1, count.covered)
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
		xlab="Number of Samples", ylab="Number of Genes Covered",
		type="l", col="blue", ylim=c(0,max(length(total.genes))),
		main=paste("Commonly Covered Genes (total = ",length(total.genes),")",sep=""))
abline(h=length(total.genes),col="red")
dev.off()