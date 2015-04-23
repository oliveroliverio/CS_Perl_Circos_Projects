#!/usr/bin/perl

while(<>) {
	@a = split;
	($mir, $exp) = ($a[0], $a[1]);
	$mir =~ s/MIR/miR-/;
	print "$mir\t$exp\n";
}

