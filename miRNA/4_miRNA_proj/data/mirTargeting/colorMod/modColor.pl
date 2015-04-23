#!/usr/bin/perl
use English; use strict; use warnings;
my(@a, @b);
my($line1, $line2, $line3, $type);

while(<DATA>) {
	$line1 = $_; # first line
	chomp($line1);
	$line2 = <DATA>; # second line  ::: apply this 
					# for any number of lines you want
	chomp($line2);
	$line3 = <DATA>;
	chomp($line3);
	@a = split(/\t/, $line1);
	($type) = ($a[4]);
	@b = split(/\t/, $line2);
	print "$line1\n";
	print "$line2\t$type\n";  # paste type to the 2nd line.  
	print "$line3\t$type\n";
	####### WORKS!!!!! #############	
}


__DATA__
116	miRNA	42	43	type=4
116	Genes	129	130	
116	stuff	3	4
11	miRNA	9	10	type=7
11	Genes	168	169
11	stuff 	4 	5
8	miRNA	9	10	type=7
8	Genes	138	139
8 	stuff	6	8
