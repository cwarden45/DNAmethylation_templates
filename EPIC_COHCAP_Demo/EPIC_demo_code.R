#setwd("C:\\Users\\Charles\\Documents\\COHCAP\\DNAmethylation_templates\\EPIC_COHCAP_Demo")
library(minfi)
library(COHCAP)

#read using minfi
description.file = "minfi_sample_description.txt"

meta.table = read.table(description.file, head=F, sep="\t")
RG.raw = read.metharray(basenames=paste("Demo\ Data\ EPIC",meta.table[,1],sep="/"))
methyl.norm = preprocessIllumina(RG.raw, bg.correct = TRUE, normalize = "controls", reference = 1)
beta.table = getBeta(methyl.norm)
probes = rownames(beta.table)

output.table = data.frame(SiteID=probes, beta.table)
beta.file = "minfi.txt"
write.table(output.table, file=beta.file, sep="\t", quote=F, row.names=F)

#analyze using COHCAP
island.file = "additional_annotations/COHCAP_EPIC_UCSC_TSS1500_plus_1st_Exon.txt"
sample.file <- "COHCAP_sample_description.txt"
project.folder <- getwd()
project.name <- "EPIC_Test"

beta.table <- COHCAP.annotate(beta.file, project.name, project.folder, platform="custom", annotation.file = island.file)
COHCAP.qc(sample.file, beta.table, project.name, project.folder)

#one-vs-one comparison
sample.file <- "COHCAP_paired_description.txt"