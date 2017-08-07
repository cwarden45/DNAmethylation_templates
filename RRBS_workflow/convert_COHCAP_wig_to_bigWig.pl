use warnings;
use strict;

my $wig_folder = "../Result/COHCAP_Results/CpG_Site/[compID]_wig";
my $wig_binary = "/opt/wigToBigWig";
my $genome = "hg19";

opendir DH, $wig_folder or die "Failed to open $wig_folder: $!";
my @files= readdir(DH);

foreach my $file (@files){
	if($file =~ /.wig$/){
		my $wig = "$wig_folder/$file";
		my $bigWig = $wig;
		$bigWig =~ s/.wig$/.bigwig/;
		if(!(-f($bigWig))){
			print "Converting $wig...\n";
			my $command = "$wig_binary $wig http://hgdownload.cse.ucsc.edu/goldenPath/$genome/bigZips/$genome.chrom.sizes $bigWig";
			system($command);
		}#end if(!(exists($bigWig)))
	}#end if($file =~ /.wig$/)
}#end foreach my $file (@files)

closedir(DH);
exit;