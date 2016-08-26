### Order to Run Scripts ###

1) minfi_preprocess.R


### Dependencies (some optional) ###

minfi: http://bioconductor.org/packages/release/bioc/html/minfi.html

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name | Name of differential expression comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to plot multiple groups|
|dmr_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|sample_description_file|Name of Sample Description File|
|beta_prefix|Table of beta values will be called *[beta_prefix]_[beta.normalization].txt*|
|beta_normalization|Normalization used to derive beta values.  Can be "illumina", "funnorm", or "noob"|
