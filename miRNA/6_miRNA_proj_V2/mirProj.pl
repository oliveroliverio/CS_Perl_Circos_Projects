#!/usr/bin/perl
# from beginning

use English;
use warnings;
use Data::Table;
use strict;

#######--------------Part 1.----------#########
# Read in the giant database and crosserefence DE genes
# to these in file. 
# unless exists the formatted file from the big one. 
my $bigDB = $ARGV[0];
open(IN0, $bigDB);



##################-----------Open Probe-Gene file  ------------------
# Read in the probe_gene2.txt file, which contains the probeID matched
# with Gene.  Put it in pgHash, which we'll use to map back to the 
# 3cellytpes file.  This file has expression values for the 3 cell types
# but they only have the probe IDs.  We need to substitute those IDs with
# gene symbols. 
our %pgHash = ();
my $probe_gene = $ARGV[1];
open(IN3, $probe_gene);
while(<IN3>) {
	next if $_ =~ /#/;
	my @a = split;
	(my $probeID, my $geneID) = ($a[0], $a[1]);
	$geneID =~ s/(\/\/.+)//g;
	$pgHash{$probeID} = $geneID;
}
close IN3;

# foreach my $i (keys %pgHash) {
# 	print "$i\t$pgHash{$i}\n";
# }


################   BACKBURNER #### open 3 files instead
##########-----------Open DE genes file (all 3 cell types)
# Substitute probe ID with Gene symbol and append cell type.
system "rm DE_expression_3_celltypes.txt";
my $de_genes = $ARGV[2];
open (OUT1, ">DE_expression_3_celltypes.txt");
open(IN1, $de_genes);  
{
	my $data = undef;
	my $header = undef;
	local $/ = ">";
	while(<IN1>) {
		next if /#/; # skips the comments in the file
		/\n/;
		$header = $PREMATCH;
		$data = $POSTMATCH;
		my ($celltype, $m) = split /\s/, $header;
		my @s = split(/\n/, $data); 
		# $data is a huge chunk, and since the input record
		# separator is no longer \n, you have to use \n as 
		# the delimiter.
		while ($header =~ /\S+/g) { 
			# had a massive problem with header ">" being printed twice
			# so to avoid this, make condition such that header must have
			# something in it, then print.  
			print OUT1 ">$header\n";  
			for my $i (0..$#s) {
				my $line = $s[$i];
				my ($id, $exp) = split /\s/, $line;
				my $subbed = &subProbeID($id);
				next if ($subbed eq "X");
				print OUT1 "$subbed\t$exp\n";
			}
		}
	}
}
close OUT1;



###########------------Open DEmirs file
# my $de_mirs = $ARGV[3];
# open(IN2, $de_mirs);
# while(<IN2>, $de_mirs) {
# 	print "$_\n";
# }





#---------------Subroutines------------- 	
#-------- Used in "Open probe gene file" section-----  
sub subProbeID () {
	my $id = shift;
	if (exists $pgHash{$id}) {$pgHash{$id}}
	else {"X"}

}


#---------------- System commands-------

#system "more delete.txt";









