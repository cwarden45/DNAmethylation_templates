use warnings;
use strict;

my $max_job_count = 10;
my %finished_samples=();

my $parameter_file = "parameters.txt";
my $ref = "";
my $alignment_folder = "";
my $reads_folder = "";
my $threads="";
my $email="";
my $strand_type = "";
my $trim_galore_path = "";
my $Bismark_path = "";
my $min_length = "";
my $cov_cutoff = "";
my $quant_method = "";
my $job_count = 0;

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
		
		if($param eq "Bismark_Ref"){
			$ref=$value;
		}elsif($param eq "Alignment_Folder"){
			$alignment_folder=$value;
		}elsif($param eq "Reads_Folder"){
			$reads_folder=$value;
		}elsif($param eq "Threads"){
			$threads=$value;
		}elsif($param eq "Cluster_Email"){
			$email=$value;
		}elsif($param eq "Strand"){
			$strand_type=$value;
		}elsif($param eq "Trim_Galore_Path"){
			$trim_galore_path=$value;
		}elsif($param eq "Bismark_Path"){
			$Bismark_path=$value;
		}elsif($param eq "Min_Trim_Length"){
			$min_length=$value;
		}elsif($param eq "Min_Coverage"){
			$cov_cutoff=$value;
		}elsif($param eq "Quantification_Method"){
			$quant_method=$value;
		}
		
	}#end if($line_count > 1)
}#end while(<PARAM>)

close(PARAM);

if(($quant_method eq "")||($quant_method eq "[[required]]")){
	die("Need to enter a value for 'Quantification_Method'!\n");
}

if(($ref eq "")||($ref eq "[[required]]")){
	die("Need to enter a value for 'Bismark_Ref'!\n");
}

if(($alignment_folder eq "")||($alignment_folder eq "[[required]]")){
	die("Need to enter a value for 'Alignment_Folder'!\n");
}

if(($reads_folder eq "")||($reads_folder eq "[[required]]")){
	die("Need to enter a value for 'Reads_Folder'!\n");
}

if(($threads eq "")||($threads eq "[[required]]")){
	die("Need to enter a value for 'Threads'!\n");
}

if(($email eq "")||($email eq "[[required]]")){
	die("Need to enter a value for 'Cluster_Email'!\n");
}

if(($strand_type eq "")||($strand_type eq "[[required]]")){
	die("Need to enter a value for 'Strand'!\n");
}

if(($trim_galore_path eq "")||($trim_galore_path eq "[[required]]")){
	die("Need to enter a value for 'Trim_Galore_Path'!\n");
}

if(($Bismark_path eq "")||($Bismark_path eq "[[required]]")){
	die("Need to enter a value for 'Bismark_Path'!\n");
}

if(($min_length eq "")||($min_length eq "[[required]]")){
	die("Need to enter a value for 'Min_Trim_Length'!\n");
}

if(($cov_cutoff eq "")||($cov_cutoff eq "[[required]]")){
	die("Need to enter a value for 'Min_Coverage'!\n");
}

my $trimmed_folder = "$reads_folder/trimmed_reads";
my $command = "mkdir $trimmed_folder";
system($command);

opendir DH, $reads_folder or die "Failed to open $reads_folder: $!";
my @files= readdir(DH);

foreach my $file (@files){
	if($file =~ /(.*)_S\d+_L\d{3}_R1_001.fastq$/){
		my ($sample) = ($file =~ /(.*)_S\d+_L\d{3}_R1_001.fastq$/);
		my $userBam = "$alignment_folder/$sample.bam";
		if((!exists($finished_samples{$sample}))&($job_count < $max_job_count)&!(-f $userBam)){
			print "$sample\n";
			
			$job_count++;
			
			my $read1 = "$reads_folder/$file";
			my $read2 = $read1;
			$read2 =~ s/_R1_001.fastq/_R2_001.fastq/g;
			
			my $sample_shell = "$sample\_preprocess.sh";
			open(OUT, "> $sample_shell") || die("Could not open $sample_shell!");

			print OUT "#!/bin/bash\n";
			print OUT "#\$ -M $email\n";
			print OUT "#\$ -m bea\n";
			print OUT "#\$ -N bis$job_count\n";
			print OUT "#\$ -q all.q\n";
			print OUT "#\$ -pe shared $threads\n";
			print OUT "#\$ -l vf=10G\n";
			print OUT "#\$ -j yes\n";
			print OUT "#\$ -o bis$job_count.log\n";
			print OUT "#\$ -cwd\n";
			print OUT "#\$ -V\n";	
			
			#trim reads
			my ($tg_prefix) = ($file =~ /(.*).fastq$/);
			my $trim1 = "$trimmed_folder/$tg_prefix\_val_1.fq";
			my $trim2 = $trim1;
			$trim2 =~ s/_R1_001_val_1.fq/_R2_001_val_2.fq/g;
			print OUT "$trim_galore_path/trim_galore --paired --rrbs $read1 $read2 -o $trimmed_folder --length $min_length\n";
			
			#align reads
			my $outputfolder = "$alignment_folder/$sample";
			print OUT "mkdir $outputfolder\n";
			
			if ($strand_type eq "none"){
				print OUT "$Bismark_path/bismark --non_directional -o $outputfolder --genome_folder $ref -p $threads -1 $trim1 -2 $trim2\n";
			}elsif($strand_type eq "yes"){
				print OUT "$Bismark_path/bismark -o $outputfolder --genome_folder $ref -p $threads -1 $trim1 -2 $trim2\n";
			}elsif($strand_type eq "reverse"){
				print OUT "$Bismark_path/bismark --pbat -o $outputfolder --genome_folder $ref -p $threads -1 $trim1 -2 $trim2\n";
			}else{
				print "Make rule for strand $strand_type\n";
				exit;
			}
			#methylation calls and summary
			my $bismarkBam = "$outputfolder/$tg_prefix\_val_1_bismark_bt2_pe.bam";
			
			if($quant_method eq "Bismark"){
				#methylation calls and summary
				#counts are in .bismark.cov file
				print OUT "$Bismark_path/bismark_methylation_extractor -o $outputfolder --buffer_size 10G --comprehensive --bedGraph --counts --no_overlap --report --genome_folder $ref -p $bismarkBam\n";

				#remove other files
				print OUT "rm $outputfolder/CHG_context_$tg_prefix\_val_1_bismark_bt2_pe.txt\n";
				print OUT "rm $outputfolder/CHH_context_$tg_prefix\_val_1_bismark_bt2_pe.txt\n";
				print OUT "rm $outputfolder/CpG_context_$tg_prefix\_val_1_bismark_bt2_pe.txt\n";				

				#get strand information (but includes all sites, even without coverage)
				my $cov_file = $bismarkBam;
				$cov_file =~ s/.bam$/.bismark.cov.gz/;
				my $stranded_file = $cov_file;
				$stranded_file =~ s/.cov.gz/.cov.v2.txt/;
				#print OUT "$Bismark_path/coverage2cytosine -o $stranded_file --genome_folder $ref $cov_file\n";
			}

			#sort and index bam for visualization and/or methylKit quantification
			my $sortPrefix = $userBam;
			$sortPrefix =~ s/.bam$//g;
			print OUT "samtools sort $bismarkBam $sortPrefix\n";
			print OUT "samtools index $userBam\n";
			
			#compress reads
			print OUT "gzip $read1\n";
			print OUT "gzip $read2\n";

			print OUT "gzip $trim1\n";
			print OUT "gzip $trim2\n";
			
			#submit job
			$command = "qsub $sample_shell";
			system($command);
		}#end if(!exists($finished_hash{$sample}))
	}#end if($file =~ /(.*)_S\d+_L\d{3}_R1_001.fastq/)
}#end foreach my $file (@files)
exit;
