# extractseq.pl
# This script will extract the fasta sequence and input to database 
# mysql -u k_luo -h mysql-sandbox.igsp.duke.edu -D kl_yeast -p
# Enter password: Ysctp5

#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;

my ($seqfile, $chrNum);
my @chrCha;
GetOptions("seqfile=s" => \$seqfile, "chrNum=i" => \$chrNum);
open(FILE, $seqfile) or die "unable to open the file";
my $directory = "/home/kl124/MNaseData/Model/PWM/genomeseq/output/";
my $outputfile_seq = $directory."sequence"."Chr".$chrNum;
my $outputfile_chr = $directory."chr".$chrNum;

open(OUT_seq, ">$outputfile_seq") or die "Cannot open file \"$outputfile_seq\" to write to!\n";
open(OUT_chr, ">$outputfile_chr") or die "Cannot open file \"$outputfile_chr\" to write to!\n";


my $seq = '';
my $chrom = '';
my $count=0;

while (<FILE>) {
   	my $line = $_;
	# discard blank line
	if ($line =~ /^\s*$/) {
		next;
	# discard comment line
	} elsif($line =~ /^\s*#/) {
		next;
	# read fasta header line to chrom
	} elsif($line =~ /^>/) {
		$chrom = $line;
		$chrom =~ s/>//g;
		next;
	# keep line, add to sequence string
	} else {
		$seq .= $line;
	}
}

# remove non-sequence data (whitespace)from $sequence string
$seq =~ s/\s//g;

print OUT_seq $seq;
print OUT_chr $chrom;
$count = length($seq);
print "\nNumber of base pairs in chr ". $chrNum. " is ". $count."\n"; 
 
@chrCha = split(/chr/, $chrom);
print @chrCha;

my $DataBaseName    =   "kl_yeast";
my $DataBaseHost    =   "mysql-sandbox.igsp.duke.edu";
my $DataBaseUser    =   "k_luo";
my $DataBasePass    =   "Ysctp5";

#$dbh = Database Handle
my $dbh =   DBI->connect("DBI:mysql:database=$DataBaseName;host=$DataBaseHost",
                                     "$DataBaseUser", 
                                     "$DataBasePass",
                             	 { RaiseError => 1, AutoCommit => 0 }) 
			or die "Unable to connect to $DataBaseHost because $DBI::errstr";

# Read the mysql_serverinfo and mysql_stat database handle attributes
my $ServerInfo  =  $dbh->{'mysql_serverinfo'};
print "Server Info: $ServerInfo\n";

# Create the table using the SQL CREATE command, with the above structure
$dbh->do("CREATE TABLE IF NOT EXISTS sacCer2_genome(chrom tinyint(3) not null, seq longtext not null)");

# Write a record to the table                                       
$dbh->do("INSERT INTO sacCer2_genome (chrom,seq) VALUES ('$chrNum', '$seq')");

$dbh->disconnect();


close FILE;
close OUT_chr;
close OUT_seq;

exit;

