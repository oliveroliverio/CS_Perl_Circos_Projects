while(<>) {
	chomp;
	$line = $_;
	$line =~ s/\/\/.+//;
	print "$line\n";
}
