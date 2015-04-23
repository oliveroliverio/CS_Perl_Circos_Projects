#!/usr/bin/perl
# This program reads in the miRNA locations file miR_loc.txt
# and has output output.txt
# Usage:  perl foo.pl


open (OUT, ">textID_mir.txt");
while(<>) {
	next if $_ =~ /#/; # skips the header stuff in the file
	@d = split;
	($chr, $start, $stop, $name) = ($d[0], $d[3], $d[4], $d[8]);
	@a = split(';', $name);
	($mir) = $a[2] =~ /Name=(.+)/;
	($mir) = $mir =~ /hsa-(.+)/; # get rid of 'hsa-'
	$chr =~ s/chr/hs/; # convert chr names to hs for circos
	print OUT "$chr\t$start\t$stop\t$mir\tid=2\n";
	
}
system "more textID_mir.txt";
