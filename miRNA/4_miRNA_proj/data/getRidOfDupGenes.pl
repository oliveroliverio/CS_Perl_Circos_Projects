#!/usr/bin/perl
# now that I have a list of mir targets of DE genes, I need to 
# make a list of genes label plot, so I need to get 
# get rid of duplicates
# usage:  perl getRidDupGenes.pl mir_gene_targetlist.txt > genesNoDup.txt

while(<>) {
	chomp;
	@a = split;
	$gene = $a[1];
	push(@b, $gene);
}

foreach $i (@b) {
	next if $seen{$i}++;
	push @unique, $i;
}

for $q (0..$#unique) {
	print "$unique[$q]\n";
}
