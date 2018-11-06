### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.

### Order to Run Scripts ###

1) `minfi_preprocess.R` or `reformat_GenomeStudio_Pval_01.R`

2) `qc.R`

3) `run_COHCAP.R`

Most results will be provided in the current working directory of the script.

### Dependencies (some optional) ###

**minfi**: http://bioconductor.org/packages/release/bioc/html/minfi.html

**BMIQ Normalization**:https://code.google.com/archive/p/bmiq/

**Island/Promoter Mapping Files**: https://sourceforge.net/projects/cohcap/files/additional_Bioconductor_annotations.zip
- Please note that "promoter" includes probes downstream of TSS (in 1st exon and 5'UTR)

**COHCAP**: https://bioconductor.org/packages/release/bioc/html/COHCAP.html

**RColorBrewer**: https://cran.r-project.org/web/packages/RColorBrewer/

### Parameter Values ###
| Parameter | Value|
|---|---|
|comp_name | Name of differential methylation comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to plot multiple groups|
|plot_type | Are are QC and heatmap labels "discrete" or "continous"?  Use commas to describe multiple variables.  If continous, orange=high, green=low|
|dmr_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|sample_description_file|Name of Sample Description File|
|project_folder|Path to output folder for selected, final results|
|beta_prefix|Table of beta values will be called *[beta_prefix]_[beta.normalization].txt*|
|beta_normalization|Normalization used to derive beta values.  Can be "illumina", "funnorm", or "noob"|
|island_mapping|Maps probeIDs to location, gene, and island|
|promoter_mapping|Maps probeIDs to location, gene, and promoter|
|cluster_distance| Distance metric for QC dendrogram.  Can be *Euclidean* or *Pearson_Dissimilarity*|
|treatment_group|Treatment group for primary variable; enter *continuous* for a continuous variable and a correlation will be provided|
|site_pvalue|Maximum p-value for filtering CpG sites|
|site_fdr|Maximum p-value for filtering CpG sites|
|site_delta_beta|Absolute value of minimum delta-beta value for filtering CpG sites|
|methyl_cutoff|Threshold for identifying a site/island as methylated|
|unmethyl_cutoff|Threshold for identifying a site/island as unmethylated|
|island_pvalue|Maximum p-value for filtering CpG sites|
|island_fdr|Maximum p-value for filtering CpG sites|
|island_delta_beta|Absolute value of minimum delta-beta value for filtering CpG sites|
|COHCAP_output_format|Either "xls" or "txt"|
|COHCAP_num_groups|Number of groups for COHCAP workflow|
|max_cluster_dist|If refining annotations in COHCAP.avg.by.island(), this is the maximum distance between adjacent filtered sites|
|min_sites_per_island|Minimum number of probes to for differentially methylated island/promoter|
|COHCAP_paired|Should COHCAP look for a second variable to consider sample pairing?|
|wig_types|Can be set to "avg" (group level and delta-beta, for 2-groups), "sample" (separate .wig file per sample", "avg.and.sample", or "none"
|expression_file|If integrating with gene expression, table of values formatted for appropriate workflow.  Otherwise, "NULL"|
