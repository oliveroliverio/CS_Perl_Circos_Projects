#!/usr/bin/perl
# This takes a regular list and turns it into format
# circos understands

while(<>) {
	chomp;
	$mir = $_;
	push(@a, $mir);
}

for $i (0..$#a) {
	print "miRNA\t$i\t$i\t0\tname=$a[$i]\ttype=1\n";
}
#$length =  $#a;
#print $length;

#for $i (0..$#a) {
#	print "$a[$i]\n";
#}
