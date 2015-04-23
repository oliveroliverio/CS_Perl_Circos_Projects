#!/usr/bin/perl
# from beginning

use English;
use warnings;
use Data::Table;
use strict;
use Data::Printer;
#use Syntax::Collector;
#use List:MoreUtils;
use Util::Any;

# usage (Testing)
# perl projMir.pl smallerDB.txt 2_<tab> 3_<tab>

# usage (For Realz)
# perl projMir.pl 1_<tab> 2_<tab> 3_<tab> 


system "clear";
#############################------GENE SIDE-----###############

##################-----------Open Probe-Gene file  ------------------
# Read in the probe_gene2.txt file, which contains the probeID matched
# with Gene.  Put it in pgHash, which we'll use to map back to the 
# 3cellytpes file.  This file has expression values for the 3 cell types
# but they only have the probe IDs.  We need to substitute those IDs with
# gene symbols. 
our %pgHash = ();
my $probe_gene = "<2_probeGeneMap.txt";
open(IN3, $probe_gene);
while(<IN3>) {
	next if $_ =~ /#/;
	my @a = split;
	(my $probeID, my $geneID) = ($a[0], $a[1]);
	$geneID =~ s/(\/\/.+)//g;
	$pgHash{$probeID} = $geneID;
}
close IN3;

##########-------<START>---Open DE genes file (all 3 cell types)
# Substitute probe ID with Gene symbol and append cell type.
system "rm DE_expression_3_celltypes.txt";
my $de_genes = "<3_DE_genes.txt";
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
		my ($celltype, $m) = split(/\s/, $header);
		my @s = split(/\n/, $data); 
		# $data is a huge chunk, and since the input record
		# separator is no longer \n, you have to use \n as 
		# the delimiter.
		if ($header =~ /\S+/g) { 
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
close IN1;

############ <END> DE gene formating done ##########

########### <START> Mapping genes to the mirs targeted to them  ######
#	I can't make a hash from this, each key must 
#	must be unique.  In this dataset, one mir can target many genes
# 	I need to construct a has of an array.  
my $bigDB = "<smallerDB.txt";
open(IN0, $bigDB);
our %mir_gene_hashArray = ();
while(<IN0>) {
	my $line = $_;
	next if $line =~ /#/;
	my @a = split;
	my($mir, $symbol) = ($a[1], $a[3]);
	($mir) = $mir =~ /hsa-(mi[Rr]-\d+)/; # get rid of hsa and other shit
	# $mir_gene_hash{$mir} = $symbol; Can't do this, I'll be missing connections
	push(@{$mir_gene_hashArray{$mir}}, $symbol);
}


# my $gene = undef;
# print "Hash content\n";
# foreach my $k (keys %mir_gene_hashArray) {
# 	print "$k\n";
# 	foreach (@{$mir_gene_hashArray{$k}}) {
# 		$gene = $_; # $_ is the contents in each key of hash $k
# 		print "$gene\t";
# 	}
# 	print "\n\n";
# }


######### <START> Map above genes to our expressed genes ############
###				This will make our DE gene track	   ############
###				ExpMirs::DE_genes					  #####
my $deGenesFile = "<DE_expression_3_celltypes.txt";
open(OUT2, ">deGenes_targeted.txt");
open(IN2, $deGenesFile);
{
	my($header1, $data1, $celltype, $m) = undef;
	local $/ = ">";
	while(<IN2>) {
		next if /#/;
		/\n/;
		$header1 = $PREMATCH;
		$data1 = $POSTMATCH;
		my ($celltype, $m) = split(/\t/, $header1);
		my @s = split(/\n/, $data1); 
		if ($header1 =~ /\S+/g) { 
			print OUT2 "> $header1\n";
			for my $i (0..$#s) {
				my ($gene, $expr) = split(/\t/, $s[$i]);
				## throw this gene into subfunction, see
				## if exists in hash of arrays
				my $mir = &checkGene_in_expMir($gene);
				print "$mir\n";
			}
		}		
	}
}
close OUT2;
close IN2;

###########------------Open DEmirs file
# my $de_mirs = $ARGV[3];
# open(IN2, $de_mirs);
# while(<IN2>, $de_mirs) {
# 	print "$_\n";
# }




#---------------Subroutines------------- 	
#-------- Used in "Open probe gene file" section-----  
sub subProbeID() {
	my $id = shift;
	if (exists $pgHash{$id}) {$pgHash{$id}}
	else {"X"}

}

sub checkGene_in_expMir() {
	my $gene = shift;
	print "Is gene $gene in here?\n";
	foreach my $k (keys %mir_gene_hashArray) {
		print "$k\n";
		foreach(@{$mir_gene_hashArray{$k}}) {
			my $gene2 = $_;
			print "Gene in hasharray: $gene2\n";
			if($gene2 eq $gene) {
				return $k;
				last;
			}
			else {
				next;
			}
		}
	}
	print "Loop ended, no gene found\n";
}


#---------------- System commands-------

#system "more delete.txt";



#2do 
#  look up one to many relationship in perl 

# make a table like below:
# > Expressed Gene track:  iCTB
# mir	mir_loc		mir_exp 	gene	gene_loc	gene_exp


############## useful code snippets ################
###################################################
#####<See if value in hash exists>#################

#     I can't do the original hash{key} = value sub
#	  Like from the translate script.  So this is what 
#	  I've come up with

## 1. Put values of hash into array.  This is useful for mapping
## 	  genes in our DE list to genes in our mirTarget list

# my @a = values %mir_gene_hash;
# for my $i (0..$#a) {
# 	print "$a[$i]\n";
# }

#  2. Turn that array back into a hash, do exists test;
#	
# 	my %a = map {$_ => 1} @a;
#	if(exists$a{$queryGene}) {...}

# example: our @mirGeneHashArray = values %mir_gene_hash;
#		   our %mirGeneHashValues = map {$_ => 1} @mirGeneHashArray;
#



########### processing line by line ###########
# use warnings; use strict; use English;

# system "clear";
# my($id, $miRNA, $type, $line1, $line2) = undef;
# my(@a,@b);
# while(<DATA>) {
# 	chomp;
# 	$line1 = $_;
# 	chomp($line1);
# 	$line2 = <DATA>;
# 	chomp($line2);
# 	@a = split(/\t/, $line1);		
# 	($type) = ($a[4]);
# 	print "$line1\n";
# 	print "$line2\t$type\n";

# }

# __DATA__
# 116	miRNA	42	43	type=4
# 116	Genes	129	130
# 11	miRNA	9	10	type=7
# 11	Genes	168	169



############## has of arrays ##########
# my $gene = undef;
# print "Hash content\n";
# foreach my $k (keys %mir_gene_hashArray) {
# 	print "$k\n";
# 	foreach (@{$mir_gene_hashArray{$k}}) {
# 		$gene = $_; # $_ is the contents in each key of hash $k
# 		print "$gene\t";
# 	}
# 	print "\n\n";
# }



