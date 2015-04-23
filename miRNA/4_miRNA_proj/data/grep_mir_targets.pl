#!/usr/bin/perl
# Script takes in a list of miR names and greps them against the 
# giant mir_target list to get miR targets

open(IN1, $ARGV[0]);
while(<IN1>) {
	chomp;
	my $mir = $_;
	push(@a, $mir);
}
close IN1;

system "rm output.txt";
open(OUT, ">mirs_gene_targetlist.txt");
open(IN2, $ARGV[1]);
while(<IN2>) {
	next if $_ =~ /#/; # skip header stuff.  
	chomp;
	$line = $_;
	for $i (0..$#a) {
		my $mir = $a[$i];
		if ($line =~ /$mir/) {
			@q = split(/\t/, $line);
			($mir1, $geneTarget) = ($q[1], $q[3]);
			$mir1 =~ s/hsa-//;
			print OUT "$mir1\t$geneTarget\n";
		}
	}
}
close OUT;

system "clear";


