use warnings; use strict; use English;

system "clear";
my($id, $miRNA, $type, $line1, $line2) = undef;
my(@a,@b);
while(<DATA>) {
	chomp;
	$line1 = $_;
	chomp($line1);
	$line2 = <DATA>;
	chomp($line2);
	@a = split(/\t/, $line1);		
	($type) = ($a[4]);
	print "$line1\n";
	print "$line2\t$type\n";

}

__DATA__
116	miRNA	42	43	type=4
116	Genes	129	130
11	miRNA	9	10	type=7
11	Genes	168	169
