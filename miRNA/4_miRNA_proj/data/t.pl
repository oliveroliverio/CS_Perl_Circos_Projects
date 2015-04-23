#!/usr/bin/perl
while(<>) {
	chomp;
	next if $_ =~ /#/;
#@a = split;
#	($mir, $aCTB, $iCTB) = ($a[0], $a[1], $a[2]);
	$stuff = $_;
	push(@a, $stuff);
}

foreach $i (@a) {
	next if $seen{$i}++;
	push (@b, $i);
}

for $i (0..$#b) {
	$line = $b[$i];
	$line =~ s/MIR/miR-/;
	print "$line\n";
}	




# -------------get rid of dups script------------
#while(<>) {
#	chomp;
#	$stuff = $_;
#	push(@a, $stuff);
#}
#
#foreach $i (@a) {
#	next if $seen{$i}++;
#	push(@b, $i);
#}
#for $i (0..$#b) {
#	print "$b[$i]\n";
#}
