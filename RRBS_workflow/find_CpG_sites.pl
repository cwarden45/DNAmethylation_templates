#Code written by Charles Warden (cwarden@coh.org, x60233)

use warnings;
use strict;
use Cwd 'abs_path'; 

$| =1;

my $os = $^O;
my $os_name;
if (($os eq "MacOS")||($os eq "darwin")||($os eq "linux"))
	{
		#Mac
		$os_name = "MAC";
	}#end if ($os eq "MacOS")
elsif ($os eq "MSWin32")
	{
		#PC
		$os_name = "PC";
	}#end if ($os eq "MacOS")
else
	{
		print "Need to specify folder structure for $os!\n";
		exit;
	}#end else
	
my $refFA = $ARGV[0];
my $CpG_table = $ARGV[1];

my $current_name = "";		
my $current_seq = "";

open(OUT, "> $CpG_table")||die("Cannot open $CpG_table");

open(IN,$refFA)||die("Cannot open $refFA\n");

while(<IN>){
	my $line = $_;
	chomp $line;
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	
	if($line =~ /^>/){
		if($current_name ne ""){
			print "Searching for CpGs in $current_name...\n";
			while ($current_seq =~ m/(CG)/g){
				#position is for end of string, in 1-based coordinate
				my $start = (pos $current_seq)-1;
				my $end = $start + 1;
				print OUT "$current_name\t$start\t$end\n";
			}#end while ($current_seq =~ m/(CG)/g)
		}#end if($current_name ne "")
		$current_name = $line;
		$current_name =~ s/>//;
	
		$current_seq = "";
	}else{
		$current_seq .= uc($line);
	}
}#end while(<IN>)
close(IN);

			print "Searching for CpGs in $current_name...\n";
			while ($current_seq =~ m/(CG)/g){
				#position is for end of string, in 1-based coordinate
				my $start = (pos $current_seq)-1;
				my $end = $start + 1;
				print OUT "$current_name\t$start\t$end\n";
			}#end while ($current_seq =~ m/(CG)/g)

close(OUT);

exit;