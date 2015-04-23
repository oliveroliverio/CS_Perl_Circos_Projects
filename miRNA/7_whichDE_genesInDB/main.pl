use English;  
#use diagnostics; 
use warnings; 
use Inline::Files; 
use Data::Printer;
use List::Compare;
#use strict;

my $db = "trimmedMirDB.txt";
open(DATA1, $db);
while(<DATA1>) {
	my $line = $_;
	my($mir, $gene) = split(/\t/, $line);
	($already{$gene}++) && next;
	#push(@geneArray, $gene);
	$d{$gene} = $mir;
}
close DATA1;

my $de = "DE_iCTB.txt";
open(DATA2, $de);
while(<DATA2>) {
	next if />/;
	$line = $_;
	my($gene, $m) = split(/\t/, $line);
	$d2{$gene} = $m;

my $lc = List::Compare->new([keys %d], [keys %d2]);
my @intersect = $lc->get_intersection;

if(!@intersect) {
	print "$_\n";
}






# ####################### This isn't working ##########
# ####################### use List::Compare instead ##########
# my $db = "trimmedMirDB.txt";
# open(DATA1, $db);
# while(<DATA1>) {
# 	my $line = $_;
# 	my($mir, $gene) = split(/\t/, $line);
# 	($already{$gene}++) && next;
# 	#push(@geneArray, $gene);
# 	$d{$gene} = $mir;
# }
# close DATA1;
# #my %d = map {$_ =>1} @geneArray;

# my $de = "DE_iCTB.txt";
# open(DATA2, $de);
# while(<DATA2>) {
# 	next if />/;
# 	$line = $_;
# 	my($gene, $m) = split(/\t/, $line);
# 	if(exists($d{$gene})) {
# 		push(@exists, $gene);
# 	}
# 	else {
# 		push(@nexists, $gene);
# 	}
# }
# close DATA2;

# my $nexistGene = "gene_not_exist.txt";
# open(OUT1, ">", $nexistGene);
# $exists = &removeDups(\@exists);
# @exists = @$exists;
# foreach(@exists) {
# 	my $gene = $_;
# 	print OUT1 "$gene\n";
# }
# close OUT1;
# p @exists;

# my $existGene = "gene_exist.txt";
# open(OUT2, ">", $existGene);
# $nexists = &removeDups(\@nexists);
# @nexists = @$nexists;
# foreach(@nexists) {
# 	my $gene = $_;
# 	print OUT2 "$gene\n";
# }
# close OUT2;
# p @nexists;

# system "wc -l gene_exist.txt";
# system "wc -l gene_not_exist.txt";




sub removeDups() {
	# fix this for proper dereferencing. 
	my $a = shift;
	my @new = @$a;

	foreach my $p (@new) {
		next if $seen{$p}++;
		push @unique, $p;
	}
	return \@unique;
}
# 2do;
# Write a script that shows the distribution of target genes in this dataset
# 		1. Counts the number of occurences of a gene and generates a histogram. 





## cool code snippets
######### See if value exists in array #######
# my %params = map { $_ => 1 } @badparams;

# if(exists($params{$someparam})) { ... }