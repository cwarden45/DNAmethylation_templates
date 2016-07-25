#setwd("C:\\Users\\Charles\\Documents\\COHCAP\\DNAmethylation_templates\\EPIC_COHCAP_Demo")
library(COHCAP)

beta.file = "minfi.txt"
#above from minfi

library(COHCAP)
island.file = "additional_annotations/COHCAP_EPIC_UCSC_TSS1500_plus_1st_Exon.txt"
sample.file <- "COHCAP_sample_description.txt"
project.folder <- getwd()
project.name <- "EPIC_Test"

beta.table <- COHCAP.annotate(beta.file, project.name, project.folder, platform="custom", annotation.file = island.file)
COHCAP.qc(sample.file, beta.table, project.name, project.folder)

#one-vs-one comparison
sample.file <- "COHCAP_paired_description.txt"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2")

project.name <- "EPIC_Test_v2"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1)
filtered.islands <- COHCAP.avg.by.site(filtered.sites, project.name, project.folder)

project.name <- "EPIC_Test_v3"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1, methyl.cutoff=0.3)
filtered.islands <- COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0.3)

project.name <- "EPIC_Test_v4"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1)
filtered.islands <- COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1, fdr.cutoff=1, num.sites=2)