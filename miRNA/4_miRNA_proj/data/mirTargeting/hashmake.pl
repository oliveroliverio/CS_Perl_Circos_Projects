#!/usr/bin/perl
# this script substitutes probeID w/ GeneID
# usage: perl hashmake.pl probe_gene2.txt 3celltypes.txt mirGeneTargets.txt mirNames.txt
use warnings;

##################-----------PART 1------------------
# Read in the probe_gene2.txt file, which contains the probeID matched
# with Gene.  Put it in pgHash, which we'll use to map back to the 
# 3cellytpes file.  This file has expression values for the 3 cell types
# but they only have the probe IDs.  We need to substitute those IDs with
# gene symbols. 
$probe_gene = $ARGV[0];
open(IN1, $probe_gene);
while(<IN1>) {
	next if $_ =~ /#/;
	@a = split;
	($probeID, $geneID) = ($a[0], $a[1]);
	$geneID =~ s/(\/\/.+)//g;
	$pgHash{$probeID} = $geneID;
}
close IN1;

##################-----------PART 2------------------
# This part will substitute the probe ID with gene symbol and push each
# line into array @b.  And print to outfile if desired.  
#
# This part will also build the gene track array after the substitution
# occurs.  I'll send that array to the getridofdups script (although, 
# I don't need to since venny said there are no genes in common).  

# Array @g contains the genes and expression values for all 3 cell types
open(OUT, ">hashout.txt");
print OUT "#cCTB\tM\tiCTB\tM\tsCTB\tM\n";
$celltypes = $ARGV[1];
open(IN2, $celltypes);
while(<IN2>) {
	chomp;
	next if $_ =~ /#/;
	@b = split;
	($cCTB, $cCTB_m, $iCTB, $iCTB_m, $sCTB, $sCTB_m ) = ($b[0], $b[1], $b[2], $b[3], $b[4], $b[5]);
	$cCTB =~ s/\d+/trans($&)/eg;
	$iCTB =~ s/\d+/trans($&)/eg;
	$sCTB =~ s/\d+/trans($&)/eg;
	push(@genetrack, $cCTB);
	push(@genetrack, $iCTB);
	push(@genetrack, $sCTB);
	#print OUT "# cCTB\tM\tiCTB\tM\tsCTB\tM\n";
	$line = "$cCTB\t$cCTB_m\t$iCTB\t$iCTB_m\t$sCTB\t$sCTB_m\n";
	#print OUT $line;
	push(@g, $line);
	#$line =~ s/(\d\d+[^\.\d+])/trans($&)/eg;
	#print OUT "$line\n";
}
close IN2;

##################-----------PART 3------------------
# Here I will create the genetrack file from the genetrack array made
# above.  WORKS!!!!  Plug this into circos. 
# I'm also creating a geneLoc hash.  This will hold the track location
# info where genes are the keys and locations (0 1, 2 3, 3 4) are the
# values.   
#
# Need to modify this to make tracks for the 3 cell types.  
open(OUT2, ">genetrack.txt");
my @uniqueArray = &getridofdups(@genetrack);
for my $k (0..$#uniqueArray) {
	my $gene = $uniqueArray[$k];
	my $f = $k+1;
	print OUT2 "Genes\t$k\t$f\t0\tname=$gene\ttype=1\n";
	my $loc = "$k\t$f";
	$geneLoc{$gene} = $loc;  
}
# hash works!
# foreach my $b (keys %geneLoc) {
# 	print "$b\t$geneLoc{$b}\n";
# }

#system "more -s genetrack.txt";

##################-----------PART 4------------------
# Here we're making a hash that contains the targeted genes as keys
# and the miRNAs that target them as values.  I need to make a hash 
# of this to add the miRNAs to array @g;  
# This portion just formats the mirGeneTargets.txt file, makes the 
# hash.  
# 
$mirGeneTargets = $ARGV[2];
open(IN3, $mirGeneTargets);
while(<IN3>) {
	chomp;
	next if $_ =~ /#/;
	my @a = split;
	(my $mir, my $gene) = ($a[0], $a[1]);
	$mir =~ s/(hsa-)//;
	$mir_gene = "$mir\t$gene";
	#push(@mir_gene_array, $mir_gene);
	$mirGeneHash{$gene} = $mir;
	#push(@geneTargetsArray, $gene);
}
close IN3;


##################-----------PART 5------------------
#Here's where the addition of miRNAs to their target genes happens, and
# we're using the mirGenesHash we built from previous code.  
# I go through array @g, which contains the 3 cell types, the DE genes
# and their expression values.  Each element in array @g is a line 
# separated by tabs.  The trans2 function will look at the gene symbol,
# and see if a miRNA target exists in the mirGenesHash.  If so, then it'll
# append the miRNA to the gene, if not, an "X" will be appended.  
for my $i (0..$#g) {
	my $line = $g[$i];
	my @b = split(/\t/, $line);
	(my $cCTB, my $cCTB_m, my $iCTB, my $iCTB_m, my $sCTB, my $sCTB_m ) = ($b[0], $b[1], $b[2], $b[3], $b[4], $b[5]);
	#substitute stuff here
	$cCTB =~ s/(.+)/trans2($&)/eg;
	$iCTB =~ s/(.+)/trans2($&)/eg;
	$sCTB =~ s/(.+)/trans2($&)/eg;

	my $line2 ="$cCTB\t$cCTB_m\t$iCTB\t$iCTB_m\t$sCTB\t$sCTB_m";
	push(@mirTarget_M_array, $line2);
}

##################-----------PART 6------------------
# Here I will create 3 histogram files for the 3 cell types using 
# array @g that I built from PART 2.  And the geneLoc hash I built from
# PART 3
open(OUT3, ">histogram_cCTB.txt");
open(OUT4, ">histogram_iCTB.txt");
open(OUT5, ">histogram_sCTB.txt");

for my $y (0..$#g) {
	my $line = $g[$y];
	my @b = split(/\t/, $line);
	(my $cCTB, my $cCTB_m, my $iCTB, my $iCTB_m, my $sCTB, my $sCTB_m ) = ($b[0], $b[1], $b[2], $b[3], $b[4], $b[5]);
	# 0 and 1: cCTB
	# 2 and 3: iCTB
	# 4 and 5: sCTB
	$cCTB =~ s/(.+)/hist_put_locs($&)/eg;
	print OUT3 "Genes\t$cCTB\t$cCTB_m\n";
	$iCTB =~ s/(.+)/hist_put_locs($&)/eg;
	print OUT4 "Genes\t$iCTB\t$iCTB_m\n";
	$sCTB =~ s/(.+)/hist_put_locs($&)/eg;
	print OUT5 "Genes\t$sCTB\t$sCTB_m\n";
	# my $line3 ="$cCTB\t$cCTB_m\t$iCTB\t$iCTB_m\t$sCTB\t$sCTB_m";
	# print "$line3";

}
close OUT3;
close OUT4;
close OUT5;
#system "more histogram_cCTB.txt";

##################-----------PART 7------------------
# Here I will make the connections using the %geneLoc hash from PART 3
# I only want to connect to these since these these are the only genes on
# our gene track.  I need to match these genes to their location (0 1, 2 3, etc)
# I need a format as below
#  1    miRNA       start location (0)           stop location (1)
#  1    Genes       start location  (2)          stop location (3)
# How to find associated miRNA, map these genes to the mirGeneHash from 
# PART 4
# First I read in the file that contains track location info of our mirs from 
# mirNames.txt and build a mirLoc hash that part 7.1 will go thorugh
$mirNames = $ARGV[3];
open(IN4, $mirNames);
while(<IN4>) {
	chomp;
	my $line4 = $_;
	my @n = split(/\t/, $line4);
	(my $start, my $stop, my $name) = ($n[1], $n[2], $n[4]);
	$name =~ s/(name=)//;
	#print "$name\n";
	($name2) = $name =~ /(miR-\d+)/;
	$mirLoc{$name} = "$start\t$stop";
	$locMir{"$start\t$stop"} = $name;
}

##############---------- PART7.1

## Push mirlocation and gene location to hash to get rid of x
open (OUT7, ">mirLocations.txt");
open(OUT6, ">connections.txt");
foreach my $c (keys %geneLoc) {
	$v++;
	foreach my $j (keys %mirGeneHash) {
		if ($j eq $c) {
			my $mir4 = $mirGeneHash{$j};
			$mir4 =~ s/-5p//;
			$mir4 =~ s/b//;
			$mir4 =~ s/a//;
			$mirLocation = &getMirLoc($mir4);
			$geneLocation2 = "$v\tGenes\t$geneLoc{$c}";
			$mirLocation2 = "$v\tmiRNA\t$mirLocation";
			$mirGeneLocHash{$mirLocation2} = $geneLocation2;
		}
	}	
}
close OUT7;
# print out the connections file. 
foreach my $e (keys %mirGeneLocHash) {
	next if $e =~ /X/;
	print OUT6 "$e\n$mirGeneLocHash{$e}\n";
}
close OUT6;


##################### Part 8 ##########################
# Add color links
#open the file made from part 7.  
# Tag these mirs:  433, 503, 424, 185, 486, 339, 421, 21 by making a hash
%mirColorHash = (
	'miR-433' => 'type=1', # red
	'miR-503' => 'type=2', # orange
	'miR-424' => 'type=3', # yellow
	'miR-185' => 'type=4', # green 
	'miR-486' => 'type=5', # blue
	'miR-339' => 'type=6', # purple
	'miR-421' => 'type=7', # pink
	'miR-21'=> 'type=8',
		);

open FILE, "<", "connections.txt" or die $!;
open (OUT8, ">connectionsColor.txt");
while(<FILE>) {
	chomp;
	my $line5 = $_;
	my @v = split(/\t/, $line5);
	(my $id, my $start, my $stop) = ($v[1], $v[2], $v[3]);
	if ($id eq "miRNA") {
		my $loc2 = "$start\t$stop";
		$mir5 = &appendColorToLinks($loc2);
		if (exists $mirColorHash{$mir5}) {
			$color = $mirColorHash{$mir5};
			$added2line = "$line5\t$color";
			print OUT8 "$added2line\n";
		}	
	} elsif ($id eq "Genes") {
		$added2line2 = "$line5\t$color";
		print OUT8 "$added2line2\n";
	}
	else {
		print OUT8 "$line5\n";
	}
}
close OUT8;

# ######################## Part 9 ######################
# # Make gene and mir track that only have links.  
# open FILE2, "<", "connections.txt" or die $!;
# open (OUT9, ">namesTrack.txt");
# while(<FILE2>) {
# 	chomp;
# 	my $line6 = $_;
# 	my @v = split(/\t/, $line6);
# 	(my $id, my $start, my $stop) = ($v[1], $v[2], $v[3]);
# 	my $loc2 = "$start\t$stop";
# 	if ($id eq "miRNA") {
# 		# do miRNA stuff here
# 		# use the mirColorHash to choose only these miRs.  

# 	}
# 	else {
# 		# do gene stuff here
# 	}
# }

##################### Part 10 ##########################
# Color links corresponding to cell type they target
#open the file made from part 7.  
# Tag these mirs:  433, 503, 424, 185, 486, 339, 421, 21 by making a hash
%cellColorHash = (
	'miR-433' => 'type=1', # red
	'miR-503' => 'type=2', # green
	'miR-424' => 'type=3', # blue
		);

open FILE, "<", "connections.txt" or die $!;
open (OUT9, ">connectionsCelltypeColor.txt");
while(<FILE>) {
	chomp;
	my $line5 = $_;
	my @v = split(/\t/, $line5);
	(my $id, my $start, my $stop) = ($v[1], $v[2], $v[3]);
	if ($id eq "miRNA") {
		my $loc2 = "$start\t$stop";
		$mir5 = &appendColorToLinks($loc2);
		if (exists $mirColorHash{$mir5}) {
			my $color = $mirColorHash{$mir5};
			$added2line = "$line5\t$color";
			print OUT9 "$added2line\n";
		}	
	}
	else {
		print OUT9 "$line5\n";
	}
}
close OUT9;


#######-------System Commands------------
# system "cp connectionsCelltypeColor.txt /Users/oliveroliverio/circos_practice/5_split_miRNA";
# system "cd /Users/oliveroliverio/circos_practice/5_split_miRNA";
# system "circos_test.sh";






#-----------------------------subroutines

sub trans() {
	if (exists $pgHash{$&}) {$pgHash{$&}}
}

sub trans2() {
	if (exists $mirGeneHash{$&}) {"$mirGeneHash{$&}::$&"}
	else {"X::$&"}
}
sub getridofdups () {
	my @b = @_;
	foreach my $p (@b) {
		next if $seen{$p}++;
		push @unique, $p;
	}
	@unique
}
sub hist_put_locs () {
	if (exists $geneLoc{$&}) {$geneLoc{$&}}
}

sub pairGeneMirs () {
	if (exists $mirGenesHash{$&}) {$mirGenesHash{$&}}
	else {"X"}
}

sub getMirLoc () {
	my $mir4 = shift;
	if (exists $mirLoc{$mir4}) {$mirLoc{$mir4}}
	else {"X"}
}

sub appendColorToLinks () {
	my $loc3 = shift;
	if (exists $locMir{$loc3}) {$locMir{$loc3}}
	else {"xx"}
}
#system "more -S geneTrack.txt";
#foreach $i (keys %pgHash) {
#	print "hello\t$i\t$pgHash{$i}\n";
#}

#2do:  venny these genes sets -> candidate genes.  
# perl venny algorithm 
# discrete math, algorithm: what's in common with all 3 datasets.  

# of the DE genes that I have, which of them are annotatedly
# targeted by mirs?  what are those mirs?




