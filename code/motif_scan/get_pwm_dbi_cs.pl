# get PWM scores using perl DBI
# mysql username and password are stored in .my.cnf under login directory (/home/home4/kevinluo/)

#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;

my ($motif_id,$score);
GetOptions("motif_id=s" => \$motif_id, "score=f" => \$score);
my $directory = "/home/home4/kevinluo/MNaseData/Model/PWM/motif_scan/PWMoutput/";
my $outputfile = $directory."motif".$motif_id."_score".$score.".txt";

open(OUT_record, ">$outputfile") or die "Cannot open file \"$outputfile\" to write to!\n";

my $pwm_results = '';

#$dbh = Database Handle
my $dsn = "DBI:mysql:kevinluo;mysql_read_default_file=$ENV{HOME}/.my.cnf";
my ($user,$password) = ("","");
my ($dbh,$sth,$query);

$dbh = DBI->connect($dsn,$user,$password, { RaiseError => 1 }) or die($DBI::errstr);
print "Connection established!\n"; 

# DEFINE A MySQL QUERY
my $sthSelect = $dbh->prepare('SELECT * FROM motif_scores WHERE motif_id = ? AND score>?');
$sthSelect->execute($motif_id, $score) or die("failed on motif_scores table");

#iterate over the set of selected rows
my $pwm_results_header = "chr \tstart \tend \tstrand \tscore \n";
my ($chr, $start, $end, $strand, $pwmscore);
print OUT_record $pwm_results_header;

while (my @record = $sthSelect->fetchrow_array())
{
	$chr = $record[1];
	if($record[3] == 1){
	$strand = '+';
	$start = $record[2];
	$end = $record[2] + length($record[5]) - 1;
	}elsif($record[3] == 0){
	$strand = '-';
	$start = $record[2] - length($record[5]) + 1;
	$end = $record[2];
	}
	$pwmscore = $record[4];

	$pwm_results = $chr."\t".$start."\t".$end."\t".$strand."\t".$pwmscore."\n";
	print OUT_record $pwm_results;
}

$dbh->disconnect();

# close OUT_record;
exit;

