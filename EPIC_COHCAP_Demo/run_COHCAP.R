beta.file = "minfi.txt"
#above from minfi

library(COHCAP)
island.file = "additional_annotations/COHCAP_EPIC_UCSC_TSS1500_plus_1st_Exon.txt"
sample.file = "COHCAP_sample_description.txt"
project.folder = getwd()
project.name = "EPIC_Test"

beta.table = COHCAP.annotate(beta.file, project.name, project.folder, platform="custom", annotation.file = island.file)
COHCAP.qc(sample.file, beta.table, project.name, project.folder)

#one-vs-one comparison
sample.file = "COHCAP_paired_description.txt"

Island.Count = c()
Site.Count = c()

filtered.sites = COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2")

project.name = "EPIC_Test_v2"
filtered.sites = COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1)
Site.Count[1] = nrow(filtered.sites)
filtered.islands = COHCAP.avg.by.site(filtered.sites, project.name, project.folder)
Island.Count[1] = nrow(filtered.islands)

project.name = "EPIC_Test_v3"
filtered.sites = COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1, methyl.cutoff=0.3)
Site.Count[2] = nrow(filtered.sites)
filtered.islands = COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0.3)
Island.Count[2] = nrow(filtered.islands)

project.name = "EPIC_Test_v4"
filtered.sites = COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1)
Site.Count[3] = nrow(filtered.sites)
Site.Count[4] = nrow(filtered.sites)
filtered.islands = COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1, fdr.cutoff=1)
Island.Count[3] = nrow(filtered.islands)

filtered.islands = COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1, fdr.cutoff=1, num.sites=2)
Island.Count[4] = nrow(filtered.islands)

comp.mat = data.frame(Island.Count, Site.Count)
print(comp.mat)
