while(<>) {
	chomp;
	$gene = $_;
	$gene =~ s/(\/\/.+)//g;
	print "$gene\n";
}
