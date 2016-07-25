#setwd("C:\\Users\\Charles\\Documents\\COHCAP\\DNAmethylation_templates\\EPIC_COHCAP_Demo")
library(minfi)

description.file = "minfi_sample_description.txt"

meta.table = read.table(description.file, head=F, sep="\t")
RG.raw = read.metharray(basenames=paste("Demo\ Data\ EPIC",meta.table[,1],sep="/"))
methyl.norm = preprocessIllumina(RG.raw, bg.correct = TRUE, normalize = "controls", reference = 1)
beta.table = getBeta(methyl.norm)
probes = rownames(beta.table)

output.table = data.frame(SiteID=probes, beta.table)
beta.file = "minfi.txt"
write.table(output.table, file=beta.file, sep="\t", quote=F, row.names=F)