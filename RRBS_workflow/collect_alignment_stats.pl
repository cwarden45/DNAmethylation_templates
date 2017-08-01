use warnings;
use strict;

my $parameter_file = "parameters.txt";
my $alignment_folder = "";
my $summary_file = "";
my $reads_folder = "";

open(PARAM, $parameter_file)||die("Cannot open $parameter_file\n");

my $line_count = 0;

while(<PARAM>){
	my $line = $_;
	chomp $line;
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	
	$line_count++;
	
	if($line_count > 1){
		my @line_info = split("\t", $line);
		my $param = $line_info[0];
		my $value = $line_info[1];
		
		if($param eq "Alignment_Folder"){
			$alignment_folder=$value;
		}elsif($param eq "aligned_stats_file"){
			$summary_file=$value;
		}elsif($param eq "Reads_Folder"){
			$reads_folder=$value;
		}
		
	}#end if($line_count > 1)
}#end while(<PARAM>)

close(PARAM);

if(($alignment_folder eq "")||($alignment_folder eq "[[required]]")){
	die("Need to enter a value for 'Alignment_Folder'!\n");
}

if(($summary_file eq "")||($summary_file eq "[[required]]")){
	die("Need to enter a value for 'aligned_stats_file'!\n");
}

if(($reads_folder eq "")||($reads_folder eq "[[required]]")){
	die("Need to enter a value for 'Reads_Folder'!\n");
}

my $trimmed_folder = "$reads_folder/trimmed_reads";

open(OUT, "> $summary_file")||die("Cannot open $summary_file\n");
print OUT "Sample\tTotal.Reads\tTrimmed.Reads\tTrimmed.Percent\tAligned.Reads\tAlignment.Rate\tCov.File\n";

opendir DH, $alignment_folder or die "Failed to open $alignment_folder: $!";
my @files= readdir(DH);

foreach my $file (@files){
	if(($file =~ /(.*).bam$/)&!($file =~/(.*).\d{4}.bam$/)){
		my ($sample) = ($file =~ /(.*).bam$/);
		print "$sample\n";
		
		my $sample_folder = "$alignment_folder/$sample";
		
		my $total = 0;
		my $trimmed = "";
		my $aligned = "";
		my $percent_trimmed = "";
		my $percent_aligned = "";
		my $cov_file = "NA";
		
		#get total reads from trimmed stats
		opendir DH2, $trimmed_folder or die "Failed to open $trimmed_folder: $!";
		my @files2= readdir(DH2);
		
		foreach my $file2 (@files2){
			if ($file2 =~ /$sample\_S\d+_L\d{3}_R1_001.fastq_trimming_report.txt/){
				open(IN, "$trimmed_folder/$file2")||die("Cannot open $trimmed_folder/$file2\n");
				while(<IN>){
					my $line = $_;
					chomp $line;
					$line =~ s/\n//g;
					$line =~ s/\r//g;
					
					if($line =~ /\s+Processed reads:\s+(\d+)/){
						($total) = ($line =~ /\s+Processed reads:\s+(\d+)/);
					}
				}#end while(<IN>)
				close(IN);
			}
		}#end foreach my $file2 (@files2)
		
		closedir(DH2);
		
		if($total == 0){
			die("Need to revise code to identify initial read count for $sample\n");
		}
		
		#get alignment rate (and infer trimmed rate from starting pairs)
		#and get Bismark coverage file (if it exists)
		my $Bismark_file="";
		opendir DH2, $sample_folder or die "Failed to open $sample_folder: $!";
		@files2= readdir(DH2);
		
		foreach my $file2 (@files2){
			if ($file2 =~ /_PE_report.txt/){
				$Bismark_file = "$sample_folder/$file2";
				open(IN, $Bismark_file)||die("Cannot open $Bismark_file\n");
				while(<IN>){
					my $line = $_;
					chomp $line;
					$line =~ s/\n//g;
					$line =~ s/\r//g;
					
					if($line =~ /Sequence pairs analysed in total:\t(\d+)/){
						($trimmed) = ($line =~ /Sequence pairs analysed in total:\t(\d+)/);
						$percent_trimmed = 100 * $trimmed/$total;
						$percent_trimmed=sprintf("%.2f", $percent_trimmed)."%";
					}
					if($line =~ /Number of paired-end alignments with a unique best hit:\t(\d+)/){
						($aligned) = ($line =~ /Number of paired-end alignments with a unique best hit:\t(\d+)/);
					}
					if($line =~ /Mapping efficiency:\t(\d+.\d+\%)/){
						($percent_aligned) = ($line =~ /Mapping efficiency:\t(\d+.\d+\%)/);
					}
				}#end while(<IN>)
				close(IN);
			}
			
			if($file2 =~ /_pe.bismark.cov/){
				$cov_file = "$sample_folder/$file2";
				$cov_file =~ s/.gz$//;
			}
		}#end foreach my $file2 (@files2)
		
		closedir(DH2);
		
		if ($Bismark_file eq ""){
			die("Need to revise code to identify Bismark alignment stat file for sample $sample\n");
		}
		
		print OUT "$sample\t$total\t$trimmed\t$percent_trimmed\t$aligned\t$percent_aligned\t$cov_file\n";
	}#end if($file =~ /(.*).bam$/)
}#end foreach my $file (@files)

closedir(DH);

exit;