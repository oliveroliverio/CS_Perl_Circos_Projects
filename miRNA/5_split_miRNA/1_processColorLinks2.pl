use warnings; use English; use strict;
system "clear";
my($line1,$line2) = undef;
while(<DATA>) {
	chomp;
	$line1 = <DATA>;
	$line2 = <DATA>;
	print "$line1\n";
}
__DATA__
116     miRNA   42      43      type=4
116     Genes   129     130
11      miRNA   9       10      type=7
11      Genes   168     169
8       miRNA   9       10      type=7
8       Genes   138     139


