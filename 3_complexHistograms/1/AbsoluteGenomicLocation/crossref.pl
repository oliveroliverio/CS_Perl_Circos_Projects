use English;  
#use diagnostics; 
use warnings; 
use Inline::Files; 
use Data::Printer;
use List::Compare;
#use strict;

$de = "DE_genes.txt";
open(DE, $de);
while(<DE>) {
	my $gene = $_;
	chomp($gene);
	# make DE gene hash that will be crossreferenced to big DB. 
	$hash{$gene} = $gene;
	#push(@geneArray, $gene);
}

$db = "trimmedDB_geneAbsLocation.txt";
open(DB, $db);
while(<DB>) {
	my $line = $_;
	chomp($line);
	my @a = split(/\t/, $line);
	my ($chr, $gene, $start, $stop) = ($a[0], $a[1], $a[2], $a[3]);
	#print "$gene\n";
	foreach my $i (keys %hash) {
		if ($i eq $gene) {
			print "$i\t$start\t$stop\n";
		}
	}
}
	# foreach my $i (keys %hash) {
	# 	print "$i\n";
	# 	if (exists $hash{$gene}) {
	# 		print "$i\t$start\t$stop\n";
	# 	} 
		# else {
		# 	next;
		# }
		#print "$_\n";
		# if ($_ eq "$gene") {
		# 	$hash2{$gene} = "$start\t$stop";
		# } 
		# elsif ($_ eq 'NA') {
		# 	$x++;
		# 	$hash2{"NA\t$x"} = "NA\tNA";
		# }
		# else {
		# 	next;
		# 



sub removeDups() {
	my $a = shift;
	my @new = @$a;

	foreach my $p (@new) {
		next if $seen{$p}++;
		push @unique, $p;
	}
	return @unique
}




