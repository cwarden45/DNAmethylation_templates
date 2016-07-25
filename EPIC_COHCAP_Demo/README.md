**NOTE: This script was tested with R-devel.  You might want to wait until the next Bioconductor release to test it on your own**

**1) Download and extract the Illumina EPIC Methylation Demo Dataset**

Click [here](http://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html) and select "Infinium MethylationEPIC Demo Data Set"
These are three replicates for the same sample (NA12878).

While this won't be comparible to a 2-group comparison with triplicates, it can give you can idea about how much variabilty might be expected when comparing two samples without replicates.

If COHCAP is in fact good at removing false positives, the differentially methylated region list should be small (whereas the differentially methylated site list may be longer).

**2) [Download](https://sourceforge.net/projects/cohcap/files/additional_Bioconductor_annotations.zip/download) extract EPIC custom COHCAP island annotations.**

If you wish to define your own set of annotations (in the same format), there are .bpm and .csv annotation files [here](http://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html), under "Illumina Methylation EPIC Product Files"

**3) Install [minfi](http://bioconductor.org/packages/release/bioc/html/minfi.html) and [COHCAP](https://www.bioconductor.org/packages/devel/bioc/html/COHCAP.html)**

At this point, you can run the full demo pipeline (via `Rscript minfi_R_3_3.R` followed by `Rscript COHCAP_devel_R_3_4.R`, for example), or run each command step-by-step as follows.  However, you might find that you have to use R 3.3 to get minfi to work (and you have to use R-devel for the COHCAP code, until the next Bioconductor release).

**4) Normalize the data and prepare the COHCAP input format**

'''
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
'''

**5) Annotate Probes and Produce QC Plots**

'''
library(COHCAP)
island.file = "additional_annotations/COHCAP_EPIC_UCSC_TSS1500_plus_1st_Exon.txt"
sample.file <- "COHCAP_sample_description.txt"
project.folder <- getwd()
project.name <- "EPIC_Test"

beta.table <- COHCAP.annotate(beta.file, project.name, project.folder, platform="custom", annotation.file = island.file)
COHCAP.qc(sample.file, beta.table, project.name, project.folder)
'''

**6) Filter CpG Differentially Methylated CpG Sites**

'''
sample.file <- "COHCAP_paired_description.txt"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2")
'''

*NOTE*: I have re-defined the sample.file to compare one replicate against one replicate.  Notice that I don't have to re-annotate the probes in order to run COHCAP on subsets of samples (just use a different sample description file)

Notice that there are no differentially methylated sites.  Unfortuantely, this is because we don't have replicates.  Let's try this again, removing the p-value and FDR thresholds.

'''
project.name <- "EPIC_Test_v2"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1)
'''

We can now find 70 differentially methylated sites (with no p-value or FDR filters)

**7) Identify Differentially Methylated Regions**

'''
filtered.islands <- COHCAP.avg.by.site(filtered.sites, project.name, project.folder)
'''

While we saw false positives at the site level, there are no differentially methylated regions.  Since these are repliates, that is good!

*NOTE* Here, I am using COHCAP.avg.by.site() instead of COHCAP.avg.by.island().  While COHCAP.avg.by.island() is the default stategy that works in most situations, it will not work well if you don't have replicate samples.  Importantly, with this workflow, p-values can be calculated at the region level, when no p-values could be calculated at the site level.


Just out of curiousity, let's see what happens when we set both the methylated and unmethylated thresholds to 0.3 (which is what I would recommend for clinical samples that show more heterogeniety than cell-line experiments)

'''
project.name <- "EPIC_Test_v3"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1, methyl.cutoff=0.3)
filtered.islands <- COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0.3)
'''

Now, we go from identifying 450 CpG sites to 0 differentially methylated islands.  Again, remember there should be no real differences, so defining methylation differences at the region level is less sensitive to noise than the site-level analysis.

Let's keep trying to get COHCAP to find a false positive at the differentially methylation region level, by decreasing the delta-beta threshold and getting rid of the methylated and unmethylated cutoffs

'''
project.name <- "EPIC_Test_v4"
filtered.sites <- COHCAP.site(sample.file, beta.table, project.name, project.folder, ref="rep2", fdr.cutoff=1, pvalue.cutoff=1, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1)
filtered.islands <- COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1, fdr.cutoff=1)
'''

There are still no differentially methylated islands at FDR < 0.05 (or even unadjusted p-value < 0.05 - notice the code above gets rid of the FDR requirement).  However, the number of sites per island is an important paramter (by default it is set to a minimum of four sites per island).  The code below changes that requriement to 2 sites per island.

'''
filtered.islands <- COHCAP.avg.by.site(filtered.sites, project.name, project.folder, methyl.cutoff=0, unmethyl.cutoff=1, delta.beta = 0.1, fdr.cutoff=1, num.sites=2)
'''

Now, COHCAP identifes 57 differentially methylted regions with an unadjusted p-value < 0.05.  If you leave the FDR filter at 0.05, there are still no false positives identifed at the island-level (but I wouldn't setting the minimum number of sites per island below 4).
