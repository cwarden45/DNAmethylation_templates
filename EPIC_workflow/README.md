### Order to Run Scripts ###

1) minfi_preprocess.R

2) qc.R

3a) run_COHCAP.R

or

3b) run_custom_DMR.R

Most results will be provided in the current working directory of the script.

### Dependencies (some optional) ###

**minfi**: http://bioconductor.org/packages/release/bioc/html/minfi.html

**Island/Promoter Mapping Files**: https://sourceforge.net/projects/cohcap/files/additional_Bioconductor_annotations.zip
- Strictly speaking, "promoter" can include probes downstream of TSS (such as 1st exon)

**COHCAP**: https://bioconductor.org/packages/release/bioc/html/COHCAP.html
-Currently, need to download source and install devel version locally: https://bioconductor.org/packages/devel/bioc/html/COHCAP.html

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name | Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to plot multiple groups|
|plot_type | Are are QC and heatmap labels "discrete" or "continous"?  Use commas to describe multiple variables|
|dmr_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|sample_description_file|Name of Sample Description File|
|Result_Folder|Path to output folder for selected, final results|
|beta_prefix|Table of beta values will be called *[beta_prefix]_[beta.normalization].txt*; Default = "minfi"|
|beta_normalization|Normalization used to derive beta values.  Can be "illumina", "funnorm", or "noob"; Default = "illumina"|
|island_mapping|Maps probeIDs to location, gene, and island|
|promoter_mapping|Maps probeIDs to location, gene, and promoter|
|cluster_distance| Distance metric for QC dendrogram.  Can be *Euclidean* or *Pearson_Dissimilarity*|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable and a correlation will be provided
