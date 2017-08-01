### Order to Run Scripts ###

1) `cluster_align_PE_Bismark.pl` (and, optionally, `run_methyKit.R`)

2) `collect_alignment_stats.pl`

3) `percent_methylation_table.R`

4) `qc.R`

5) `run_COHCAP.R`

Most results will be provided in the current working directory of the script.

### Dependencies (some optional) ###

**Trim Galore!**:https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

**Bismark**: https://www.bioinformatics.babraham.ac.uk/projects/bismark/

**MethylKit**: https://bioconductor.org/packages/release/bioc/html/methylKit.html

**COHCAP**: https://bioconductor.org/packages/release/bioc/html/COHCAP.html

**RColorBrewer**: https://cran.r-project.org/web/packages/RColorBrewer/

### Parameter Values ###
| Parameter | Value|
|---|---|
|Alignment_Folder|Path to Bismark Alignments|
|Reads_Folder|Path to Reads for Bismark Alignment|
|Bismark_Ref|Path to Bowtie2 Bismark Ref|
|Bismark_Path|Path to Bismark Executables|
|Trim_Galore_Path|Path to Trim Galore! Executables|
|Min_Trim_Length|Minimum Trimmed Length for Trim Galore!|
|Strand|Library type. Can be *no*, *yes*, or *reverse*.|
|Threads|Number of Threads for Bismark Alignment|
|Cluster_Email|If running alignment on a cluster, e-mail for notifications|
|Quantification_Method|Method for quantifying percent methylation / counts (can be *Bismark* or *methylKit*)|
|methyl_percent_prefix|Table of percent methylation values will be called *[methyl_percent_prefix]_[Quantification_Method].txt*|
|aligned_stats_file|Name of File to Contain Aligned Read Counts and Bismark Coverage Files (if applicable)|
|Min_Coverage|Minimum Coverage to Analyze Percent Methylation|

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
