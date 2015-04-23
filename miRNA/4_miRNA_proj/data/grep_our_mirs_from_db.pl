#!/usr/bin/perl
# this outputs a list of our DE miRNAs with the 
# formatting circos likes.  

# Usage: perl foo.pl > output.txt

while(<>) {
	chomp;
	$mir = $_;
	system "grep $mir textID_ALL_mirs.txt";
}
