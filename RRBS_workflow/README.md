### \~Temporary Note\~ ###
**I apologize for the confusion, but I would like to emphasize that these are called “templates” because I almost always have to modify the code for each project (beyond the parameter files), meaning they will be more difficult for other people to use in the same way.**  This was unfortunately not immediately clear to me when I created the templates.

I also believe that the process of writing the scripts for analysis (such as the templates) is very important for the learning process, and it is very important that you understand all the steps for analysis before presenting them in a paper.

So, I will post an update when more specific guidance / suggestions can be provided.

### Order to Run Scripts ###

1) `cluster_align_PE_Bismark.pl`

2) `collect_alignment_stats.pl`

3) `percent_methylation_table.R` or `percent_methylation_table_DESTRANDED.R`

4) `qc.R`

5) `run_COHCAP.R` or `run_methylKit.R`

Most results will be provided in the current working directory of the script.

### Dependencies (some optional) ###

**Trim Galore!**:https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

**Bismark**: https://www.bioinformatics.babraham.ac.uk/projects/bismark/

**MethylKit**: https://bioconductor.org/packages/release/bioc/html/methylKit.html

**COHCAP**: https://bioconductor.org/packages/release/bioc/html/COHCAP.html

**BiSeq**: https://bioconductor.org/packages/release/bioc/html/BiSeq.html

**RColorBrewer**: https://cran.r-project.org/web/packages/RColorBrewer/

### Parameter Values ###
| Parameter | Value|
|---|---|
|Alignment_Folder|Path to Bismark Alignments|
|Reads_Folder|Path to Reads for Bismark Alignment|
|Result_Folder|Path to output folder for selected results (with subfolders)|
|project_folder|Path for COHCAP Results|
|Bismark_Ref|Path to Bowtie2 Bismark Ref|
|Bismark_Path|Path to Bismark Executables|
|Trim_Galore_Path|Path to Trim Galore! Executables|
|Min_Trim_Length|Minimum Trimmed Length for Trim Galore!|
|Strand|Library type. Can be *no*, *yes*, or *reverse*.|
|Threads|Number of Threads for Bismark Alignment|
|Cluster_Email|If running alignment on a cluster, e-mail for notifications|
|Quantification_Method|Method for quantifying percent methylation / counts (can be *Bismark* or *methylKit*)|
|methyl_percent_prefix|Table of percent methylation values will be called *[methyl_percent_prefix]\_[Quantification_Method]\_[Min_Coverage]x.txt*|
|aligned_stats_file|Name of File to Contain Aligned Read Counts and Bismark Coverage Files (if applicable)|
|Min_Coverage|Minimum Coverage to Analyze Percent Methylation|
|Min_Coverage_Pair|If using `percent_methylation_table_DESTRANDED.R` to create a table of destranded methylation values, this is the minimum coverage for both positions in CpG.  So, if *Min_Coverage* = 10, and *Min_Coverage_Pair* = 4, then both forward and reverse position must have 4 reads and the sum of reads for both positions must be greater than 10.|
|FA_Ref|If you don't have a .bed file with CpG sites (in `percent_methylation_table_DESTRANDED.R`, this is "[genome]\_CpG.bed"), such a file can be created from a reference .fa file|
|comp_name | Name of differential methylation comparison (used to name output file)
|plot_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to plot multiple groups|
|plot_type | Are are QC and heatmap labels "discrete" or "continous"?  Use commas to describe multiple variables.  If continous, orange=high, green=low|
|dmr_groups | Names of columns in *sample_description_file* to be plotted in QC and differential methylation plots.  Use commas to include multiple variables (for multivariate model or gene list filtering)|
|sample_description_file|Name of Sample Description File|
|gene_region_BED|BED file with promoter / island locations, to create custom site mapping for covered sites.  Region names should be in the format *[gene name]\_[COHCAP region name]*.|
|gene_region_mapping|Maps probeIDs to location, gene, and promoter/island|
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
|COHCAP_paired|Should COHCAP look for a second variable to consider sample pairing?  Please note that by setting this to "TRUE", you are allowing analysis with a second variable (even if "pairing" is not 1-to-1).  Otherwise, COHCAP does not automatically recognize additional columns in the sample description file.  For some workflows, this can be set to "continuous" for analysis of a continuous co-variate.|
|wig_types|Can be set to "avg" (group level and delta-beta, for 2-groups), "sample" (separate .wig file per sample", "avg.and.sample", or "none".  Since Bismark coverage and bedGraph files already exist for RRBS samples, I recommend either "none" or "avg".|
|expression_file|If integrating with gene expression, table of values formatted for appropriate workflow.  Otherwise, "NULL"|
|genome|Name of genome build|
|Read_Pairing|Are you using single-end (*SE*) or paired-end (*PE*) reads?|
