#!/usr/bin/perl
# this script gets rid of dup miRNAs from the 
# formatted_mir_list
# usage:  foo.pl test.txt

#test.txt contains:
	# boy
    # boy
    # dog
    # cat
    # cat
    # girl 
    # apple


while(<>) {
	chomp;
	$mir = $_;
	push(@a, $mir)
}
foreach $i (@a) {
	next if $seen{$i}++;
	push @unique, $i;
}

for $q (0..$#unique) {
	print "$unique[$q]\n";
}

#my @unique = unique(@a);
#for $i (0..$#unique) {
#	print "$unique[$i]\n";
#}
#my %hash = map {$mir,1} @a;
#my @unique = keys %hash;
#
#for $i (0..$#unique) {
#	print "$unique[$i]\n";
#}
