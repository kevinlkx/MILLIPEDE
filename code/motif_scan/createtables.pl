# mysql -u k_luo -h mysql-sandbox.igsp.duke.edu -D kl_yeast -p
# Enter password: Ysctp5

#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use DBI;

my $DataBaseName    =   "kl_yeast";
my $DataBaseHost    =   "mysql-sandbox.igsp.duke.edu";
my $DataBaseUser    =   "k_luo";
my $DataBasePass    =   "Ysctp5";

#$dbh = Database Handle
my $dbh = DBI->connect("DBI:mysql:database=$DataBaseName;host=$DataBaseHost",
                                     "$DataBaseUser", 
                                     "$DataBasePass",
                             	 { RaiseError => 1, AutoCommit => 0 }) 
			or die "Unable to connect to $DataBaseHost because $DBI::errstr";

# Read the mysql_serverinfo and mysql_stat database handle attributes
my $ServerInfo  =  $dbh->{'mysql_serverinfo'};
print "Server Info: $ServerInfo\n";

# Create the table using the SQL CREATE command, with the above structure
$dbh->do("CREATE TABLE IF NOT EXISTS motif_ids(
		motif_id tinyint(3) unsigned, 
		PRIMARY KEY(motif_id),
		description text,
		background_order tinyint(4),
		minimum_score double)");

# Write a record to the table                                       

$dbh->disconnect();

exit;

