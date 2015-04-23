my $expressedMirs = "4_DE_mirs.txt";
open(MIRDATA, $expressedMirs);
while(<MIRDATA>) {
		my ($mir) = $_ =~ /(\S+[^-][\d+])/;
		print "$mir\n";
}
close MIRDATA;
