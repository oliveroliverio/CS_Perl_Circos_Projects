#!/usr/bin/perl
use warnings;
use English;

$file = "<1_hg19_predictions_S_C_aug2010.txt";
$file2 = "<asdf";
open(IN, $file);
open(OUT, ">trimmedMirDB.txt");
while(<IN>) {
	$line = $_;
	next if $line =~ /#/;
	@a = split(/\t/, $line);
	($mir, $gene) = ($a[1], $a[3]);
	$mir =~ s/(hsa-)//;
	if ($mir =~ /let/) {
		print OUT "$mir\t$gene\n";
	}
	else {
		($mir) = $mir =~ /(mi[Rr]-\d+)/;
		print OUT "$mir\t$gene\n";
	}
}
