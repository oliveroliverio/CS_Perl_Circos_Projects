#!/usr/bin/perl

while(<>) {
	@a = split;
	$mir = $a[0];
	push(@b, $mir);
}

foreach $i (@b) {
	next if $seen{$i}++;
	push(@unique, $i);
}

for $q (0..$#unique) {
	print "$unique[$q]\n";
}
