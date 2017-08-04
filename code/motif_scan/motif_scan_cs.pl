#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use DBI;
use Getopt::Long;

# Set autoflush
$| = 1;

my ($n, $motif_id, $pwm_file, $mdesc, $bfile, $min_score);

# List of options
#	n:		Background order of the Markov process
#	pwm:		PWM file in MAST format (see below)
#	motif_id:	The motif id to utilize in the mysql table
#	desc:		The description to use in the mysql table
#	bfile:		The background file of genomic percentages for different kmers
#	min_score:	The minimum motif score to enter into the database (defaults to 0)

my $result = GetOptions("n=i"=>\$n,
                        "pwm=s"=>\$pwm_file,
                        "motif_id=i"=>\$motif_id,
                        "desc=s"=>\$mdesc,
                        "bfile=s"=>\$bfile,
			"min_score:f"=>\$min_score);

if(!$motif_id || !$pwm_file || !$mdesc || !$bfile) {
  die_with_message("not all parameters present");
}
if(!defined($n)) {
  die_with_message("No n defined");
}

if(!(-e $pwm_file)) {
  die("couldn't open file $pwm_file");
}

# Establish the Perl DBI connection ($dsn, $username, $password)
######################### changed ##################################

#$dbh = Database Handle
my $dsn = "DBI:mysql:kevinluo;mysql_read_default_file=$ENV{HOME}/.my.cnf";
my ($user,$password) = ("","");
my ($dbh,$sth,$query);

$dbh = DBI->connect($dsn,$user,$password, { RaiseError => 1 }) or die($DBI::errstr);
print "Connection established!\n"; 

# my $dbh=DBI->connect("DBI:mysql:database=kl_yeast;host=mysql-sandbox.igsp.duke.edu", "k_luo", "Ysctp5");

# Format is to have two tables that hold the information
# Each motif_scan with a particular PWM is given a motif_id
# Description of this motif_id is in the motif_ids table
# The ~25 million scores per position are stored in the motif_scores table
# Table 1: motif_ids
#	Field			Type			Key
#	motif_id		tinyint(3) unsigned	PRI
#	description		text
#	background_order	tinyint(4)
#	minimum_score		double
# Table 2: motif_scores
#	Field			Type			Key
#	motif_id		tinyint(3) unsigned	PRI
#	chr			tinyint(3) unsigned	PRI
#	pos			int(11)			PRI
#	strand			tinyint(4)		PRI (1 is "+", 0 is "-")
#	score			double	
#	seq			varchar(64)			
# Table 3: sacCer2_genome
#	Field			Type
#	chrom			tinyint(3)
#	seq			longtext

# Check to make sure that the motif_id supplied is unique in the DB
######################### change ##################################
my $sth_check = $dbh->prepare("SELECT COUNT(*) FROM motif_ids WHERE motif_id=?");
my $rowval = $sth_check->execute($motif_id) || die("Failed on motif_ids table check: $DBI::errstr");
my $allrv = $sth_check->fetchall_arrayref();
if($allrv->[0]->[0] > 0) {
  die("motif_id $motif_id already exists in the DB");
}

# Get the genomic sequence for each chromosome and store it into an array called $seqs
print "Grabbing genome...\n";
######################### changed ##################################
my $sth_down = $dbh->prepare("SELECT chrom,seq FROM sacCer2_genome ORDER BY chrom");
$sth_down->execute() || die("Error selecting down sequences");
my $seqs = $sth_down->fetchall_arrayref();

# Read in the background sequence files
# Background file
#	kmer	percent
#	A	0.3090
#	C	0.1916
#	G	0.1913
#	T	0.3080
#	AA	0.3504
print "Reading in background sequences from $bfile...\n";
my %bg = ();
my $bfilein;
open($bfilein,"<$bfile") || die_with_message("Failed opening file $bfile");
# ignores any lines not of the form N+ = 0.\d+
while(my $line = <$bfilein>) {
  chomp($line);
  if($line =~ m/^([ACGT]+)\t([\d\.]+)$/) {
    # Hash the k-mer sequence to its proportion in the genome
    # "A" -> 0.3090
    # bg_all.tab has genome proportions of all k-mers from 1 to 6
    $bg{$1} = $2;
  }
}

# Lists all of the background statistics
print Dumper(\%bg);

# Read in the MAST format of the PWM
# Tab-delimited file of nucleotide percentages at each position
# ALPHABET= ACGT
# 0.3849	0.14748	0.2446	0.2230
print "reading in motif...\n";
# Read the PWM into an array
open(MIN,$pwm_file) || die("Bad motif file: $pwm_file");
my @alphabet = undef;
my %alphalup = ();
my @pos = ();
my $j = 0;
my $mwidth;
while(<MIN>) {
	chomp;
	my $placed = 0;
	my $l = $_;
	# Example: in motif_chipseq.v2.pwm header is ALPHABET= ACGT (m//i means case-insensitive)
	# Hash each nucleotide to a number (e.g. A to 1, C to 2, G to 3, T to 4)
	if($l =~ m/alphabet=\s*(\w+)/i) {
		@alphabet = split(//,$1);
		print Dumper(\@alphabet);
		my $i = 0;
		foreach my $nt (@alphabet) {
			$alphalup{$nt} = $i;
			$i++;
		}
	}
	# If the alphabet has been defined
	elsif(@alphabet) {
		my $i = 0;

		# Insert into the two-dimensional array the nucleotide percentage
		# Dimension 1: Position of the PWM
		# Dimension 2: Nucleotide (e.g. A to 1, C to 2, G to 3, T to 4)
		while($l =~ m/([0-9\-\.Inf]+)/ig) {
			$placed = 1;
			$pos[$j]->[$i] = $1;
			$i++;
		}
	}
	# All of the nucleotide percentages for position j have been read in
	# Increment j to the next nucleotide position
	if($placed) {
		$j++;
		$mwidth=$j;
	}
}

close MIN;
print Dumper(\@pos);

# Print the width of the motif
print "mwidth = $mwidth\n";
my $pos_ptr = \@pos;
my $alpha_ptr = \%alphalup;

# Enter into motif_ids a description of the motif scan
######################### change ##################################
my $sth_mid = $dbh->prepare("INSERT INTO motif_ids (motif_id,description,background_order,minimum_score) VALUES (?,?,?,?)");
$sth_mid->execute($motif_id,$mdesc,$n,$min_score) || die("Failed on motif_ids insert: $DBI::errstr");

# Prepare the statement that will enter in the motif score for each position in the genome
######################### change ##################################
my $sth_up = $dbh->prepare("INSERT INTO motif_scores (motif_id,chr,pos,strand,score,seq) VALUES (?,?,?,?,?,?)");

# Run the motif scan on each chromosome
print "Scoring sequences...\n";
my $z = 0;
# $seq->[0] has the chromosome name
# $seq->[1] has the chromosome number

# Increment over all of the chromosomes
# NOTE: pos is the beginning position of the word closest to the 5' end
# 	e.g. if the PWM is 8 bp and pos is 2300:
#		if the position is "+", then this refers to 2300-2307
#		if the position is "-", then this refers to 2293-2300

foreach my $seq (@$seqs) {
	print "chr ".$seq->[0]."\n";
	
	# Get the length of each nucleosome 
	my $len = length($seq->[1]);
	# Run the motif scan over the "+" strand  
	
	# Split the genomic sequence of the chr into an array so that each
	# position of the array holds the nucleotide sequence at that position
	my @seq = split(//,$seq->[1]);
  	
	# Increment over every position in the chromosome
	for(my $i = 0; $i < ($len - $mwidth - 1); ++$i) {
    		# Get the "word" array, which is the length of the PWM starting at
		# every position in the genome
		my @word = @seq[$i..(($i+$mwidth)-1)];
    
		# Calculate the background score
		my $bg_score = background_score(\@word,\%bg,$n);
    
		# Calculate the motif score
		my $fg_score = WM_score(\@word,$pos_ptr,$alpha_ptr);
   
		# If WM_score returns a -999, encountered an invalid nucleotide in the PWM, so skip 
		if ($fg_score == -999){
			next;
    		}
    
		# Calculate the log-odds ratio of motif score over background score
		# Since logarithm, can just substract the scores
		my $lod_score = $fg_score - $bg_score;
    
		# Only input score into the database if it meets a minimum score
		if ($lod_score > $min_score){
    			$sth_up->execute($motif_id, $seq->[0], $i+1, 1, $lod_score, join('',@word)) || 
			die("Failed on positive upload: $DBI::errstr");
    		}
  	}
  
	# Run the motif scan over the "-" strand
  	print "Reverse complementing...\n";
  	@seq = @{rev_comp(\@seq)};
  	print " ... RC complete!\n";
  	
	# Increment over every position in the chromosome
	for(my $i = 0; $i < ($len - $mwidth - 1); ++$i) {
    	
    		# Get the "word" array, which is the length of the PWM starting at
		# every position in the genome
		my @word = @seq[$i..(($i+$mwidth)-1)];

		# Calculate the background score
    		my $bg_score = background_score(\@word,\%bg,$n);
    
		# Calculate the motif score
		my $fg_score = WM_score(\@word,$pos_ptr,$alpha_ptr);
    
		# If WM_score returns a -999, encountered an invalid nucleotide in the PWM, so skip 
		if ($fg_score == -999){
			next;
    		}

		# Since logarithm, can just substract the scores
    		my $lod_score = $fg_score - $bg_score;
		
		# Only input score into the database if it meets a minimum score
	   	if ($lod_score > $min_score){
	    		$sth_up->execute($motif_id, $seq->[0], ($len-$i), 0, $lod_score, join('',@word)) 
			|| die("Failed on negative upload: $DBI::errstr");
	    	} 
  	
	}
  	
	print "Done with chr " . $seq->[0] . "\n";
}

# all done
$dbh->disconnect();


# SUBROUTINES

# Reverse-complement subroutine
sub rev_comp {
	# Assign the "+" strand sequence to word
	my $word = shift;
  	my $len = scalar(@$word);
  	my $word2 = [];
  	# Increment over each position on the chromosome
	# and find the complementary base pair
	for(my $i = 0; $i < $len; ++$i) {
    		my $nt = $word->[$i];
    		if($nt eq 'A'){
      			$word2->[(($len-1)-$i)] = 'T';
    		}
		elsif ($nt eq 'C') {
      			$word2->[(($len-1)-$i)] = 'G';
    		}
		elsif ($nt eq 'G') {
      			$word2->[(($len-1)-$i)] = 'C';
		}
		elsif ($nt eq 'T') {
      			$word2->[(($len-1)-$i)] = 'A';
    		}
		else{
      			$word2->[(($len-1)-$i)] = 'N';
    		}
  	}
  	return($word2);
}

# Motif-score subroutine
sub WM_score {
 	
	# $word is the PWM-length genomic sequence
	# $motif is the two-dimensional PWM-vector
	# $mlup is the hash that associates a nucleotide with a number
	my ($word,$motif,$mlup) = @_;
  
	my $mwidth = scalar(@$word);

	# $total is the log-probability of the given sequence matching the PWM
	my $total = 0;

	# Increment over each nucleotide in the word	
	for(my $i = 0; $i < $mwidth; ++$i) {
    		# Get the probability of a particular nucleotide at the given position	
		my $temp_nuc_score = $motif->[$i]->[$mlup->{$word->[$i]}];
    		if ($temp_nuc_score > 0){
    			$total += log($temp_nuc_score);
    		}
    		# If the PWM does not allow a nucleotide in the word at the position,
		# (if the percentage is 0 for a nucleotide), then return "-999"
		else{
			return(-999);
    		} 
  	}
  	return($total);
}

# Background-score subroutine
sub background_score {
	
	# $word is the PWM-length genomic sequence
	# $bg is the hash of a kmer sequence to its genomic percentage
	# $n is the markov order
	my ($word, $bg, $n) = @_;
  	my $mwidth = scalar(@$word);
  	my $last = 0;
  	my $total_score = 0;
  	# Increment over every position of the word
	for(; $last < $mwidth; ++$last) {
    		my $first = $last - $n;
    		# If first is less than 0, then assign it 0, otherwise, keep it at its original value
    		$first = ($first < 0) ? 0 : $first;
    		# For n = 4, starts at 0, 1, 2, 3, 4, and then increments at 4 until the end of the motif
    		# Form the subword at each position
		# 	e.g. if word is ACTTAGAT, then the scores used for a 4th order markov are
		#		genomic percentage at "A", "AC", "ACT", "ACTT", "CTTA", "TTAG", "TAGA", "AGAT"	
		my $subword = join('',@{$word}[$first..$last]);
    			$total_score += log($bg->{$subword});
  		}
  	return($total_score);
}

# Description of the script
sub die_with_message {
  my $msg = shift;
  print STDERR $msg."\n";
  print STDERR "Usage: ./motif_scan_using_nrth_order_markov_background.pl
    --motif_id=<The PK to use for motif_score table>
    --pwm=<MAST formated file describing PWM>
    --desc=<Description for motif_ids table>
    --n=<Order for Markov Background Model>
    --bfile=<Background file, as output by motif_scan_build_nth_order_markov_background.pl>
    --min_score=<Minimum log-odds score to enter into database (defaults to 0)>
    ";
  die("Execution failed");
}

