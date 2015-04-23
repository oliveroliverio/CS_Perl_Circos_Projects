#!usr/bin/perl
use warnings;

#directory: /Users/oliveroliverio/circos_practice/4_miRNA_proj/data/mirTargeting
#usage: foo.pl DE_mirs.txt hgblahblah_bigFile.txt
$file1 = $ARGV[0];

open(IN, $file1);

while(<IN>) {
	chomp;
	$mir = $_;
	push(@a, $mir); 
	# @a is an array containing all our DE miRNAs.  
}
close IN;

# for $s (0..$#a) {
# 	print "$a[$s]\n";
# }
$file2 = $ARGV[1];
open(IN2, $file2);
open(OUT, ">mirGeneTargets.txt");

while(<IN2>) { 
	# here I'm openning the massive miRNA-target
	# database, and querying it with our DE miRNAs.  
	next if $_ =~ /#/;
	chomp;
	$line = $_;
	for $i (0..$#a) {
		$mir2 = $a[$i];
		if ($line =~ /$mir2/i) { #i for case insensitive
			push(@b, $line);
			# here I'm pushing all the lines that contain
			# our DE miRNAs into @b.  
		}
	}
}
close IN2;
# below I'm getting only the miRNA names and genes 
# that it targets.  
for $x (0..$#b) {
	$dbline = $b[$x];
	my @a = split(/\t/, $dbline);
	(my $mir, my $gene) = ($a[1], $a[3]);
	print OUT "$mir\t$gene\n"; 
	push(@geneList, $gene);
	# Now I have an array conataining all my genes
	# but there may be duplicates because multiple miRNAs
	# target multiple genes.  
	
}

# This subfunction gets rid of the duplicates.
&getRidOfDups(@geneList);



#system "more -S mirGeneTargets.txt";

sub getRidOfDups { # also prints out the gene track.  
	open (OUT2, ">geneTrack.txt");
	my @b = @_;
	foreach my $i (@b) {
		next if $seen{$i}++;
		push @unique, $i;  
		# now we have an array of all the genes that are targetted
		# by our differentially expressed miRNAs.  
	}
	for my $q (0..$#unique) {
		my $z = $q + 1;
		$name = "name=" . "$unique[$q]";
		print OUT2 "Genes\t$q\t$z\t0\t$name\ttype=1\n";
	}
	close OUT2;
}

system "more geneTrack.txt";

#makes:
# gene track
# interactions list
# histogram
#-------------------

# for $x (0..$#b) {
# 	$dbLine = $b[$x];
# 	print OUT $b[$x];
# 	print OUT $dbLine;	
# 	my @a = split(/\t/, $dbLine);
# 	(my $mir, my $gene) = ($a[1], $a[3]);
# 	print OUT "$mir\t$gene\n";
# }

# system "more -S output.txt";





