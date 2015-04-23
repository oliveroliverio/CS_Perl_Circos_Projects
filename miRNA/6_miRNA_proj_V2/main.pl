use English;  
#use diagnostics; 
use warnings; 
use Inline::Files; 
use Data::Printer;
use List::Compare;
#use strict;

system "clear";
system "rm BEAST.txt";
system "rm geneNamesTrack.txt";
system "rm mirNamesTrack.txt";
system "rm connections.txt";
system "rm namesTrack2.txt";



# my($line, $mir, $gene) = undef;
my(@a);
my %mirGeneHashArray = ();
########################### -------- Open and process giant mirGene DB
my $giantDB = "trimmedMirDB.txt";
open (DATA2, $giantDB);
while(<DATA2>) {
	$line = $_;
	chomp($line);
	($mir, $gene) = split(/\t/, $line);
	#print "$mir\n";
	push(@{$mirGeneHashArray{$mir}}, $gene);
	
}
close DATA2;


#########################---------- Get only expressed mirs.  (Filter the above hash for mirs that ours express)
my $expressedMirs = "4_DE_mirs.txt";
open(MIRDATA, $expressedMirs);
while(<MIRDATA>) {
	chomp;
	my ($mir) = $_ =~ /(\S+[^-][\d+])/;
	if (exists $mirGeneHashArray{$mir}) {
		print "$mir\n";
		foreach(@{$mirGeneHashArray{$mir}}) {
			print "$_\n";
			push(@{$de_mirs_hash{$mir}}, $_);
		}
		#push(@{$de_mirs_hash{$mir}}, $mirGeneHashArray{$mir});
	}
}
close MIRDATA;

#p %de_mirs_hash;

############--------------------Parse the DE genes files (3 cell types)------

# p %mirGeneHashArray;
my $geneExpData = "DE_expression_3_celltypes.txt";
my $geneExpData2 = "DE_exp_3Celltypes_100.txt";
my $geneExpData3 = "DE_exp_3Celltypes_200.txt";
open (DATA1, $geneExpData3);
{
	local $/ = ">";
	#my($data, $header, $celltype, $m, $id ) = undef;
	my @s;
	while(<DATA1>) {
		/\n/;
		my $header = $PREMATCH;
		chomp($header);
		my $data = $POSTMATCH;
		#print "$data\n";
		my($celltype, $m) = split(/\t/, $header);
		#print "$celltype\n";
		#############-------------------------aCTB-----------
		if ($celltype =~ "aCTB") {
		# 	## create gene track for DE genes for aCTB sample
		# 	## place each cell type into respective subroutines
		# 	## use List::Compare module
			
			my @s = split(/\n/, $data);
			for my $i (0..$#s) {
				my ($id, $exp) = split(/\s/, $s[$i]);
				if ($id eq ">") {
					next;
				}
				else {
					$aCTB_hash{$id} = $exp;
				}
			}
			# pass the %aCTB_hash to subroutine that will match
			# genes in the mir-target database
			my $one_ref = &parse_aCTB(\%aCTB_hash, \%de_mirs_hash);
			
			%DE_genes_aCTB = %$one_ref;
		} 
		#############-------------------------iCTB-----------
		elsif ($celltype =~ "iCTB") {
			my @s = split(/\n/, $data);
			for my $i (0..$#s) {
				my ($id, $exp) = split(/\s/, $s[$i]);
				if ($id eq ">") {
					next;
				}
				else {
					$iCTB_hash{$id} = $exp;
				}
			}
			# pass the %aCTB_hash to subroutine that will match
			# genes in the mir-target database
			my $one_ref = &parse_iCTB(\%iCTB_hash, \%de_mirs_hash);
			
			%DE_genes_iCTB = %$one_ref;
		} 
		# ############-------------------------sCTB-----------
		else { # last cell type is sCTB
			my @s = split(/\n/, $data);
			for my $i (0..$#s) {
				my ($id, $exp) = split(/\s/, $s[$i]);
				if ($id eq ">") {
					next;
				}
				else {
					$sCTB_hash{$id} = $exp;
				}
			}
			# pass the %aCTB_hash to subroutine that will match
			# genes in the mir-target database
			my $one_ref = &parse_sCTB(\%sCTB_hash, \%de_mirs_hash);
			
			%DE_genes_sCTB = %$one_ref;
		} 
	}
}


############# Now I have all the information I need contain in this
############# hash of arrays: all genes, expression levels and 
############# their associated mirs

############# Part 2:  Make gene and mir track
# Note: these tracks are only ones w/ genes w/ connections to 
# reported mirs.  
# make sure to sort by expression value

my ($tracklength_aCTB, $aCTB_geneLoc_ref, $aCTB_hist) = &printGeneTracks(\%DE_genes_aCTB, "aCTB", 0);
%aCTB_geneLoc = %$aCTB_geneLoc_ref;
%aCTB_hist = %$aCTB_hist;

my ($tracklength_iCTB, $iCTB_geneLoc_ref, $iCTB_hist) = &printGeneTracks(\%DE_genes_iCTB, "iCTB", 0);
%iCTB_geneLoc = %$iCTB_geneLoc_ref;
%iCTB_hist = %$iCTB_hist;

my ($tracklength_sCTB, $sCTB_geneLoc_ref, $sCTB_hist) = &printGeneTracks(\%DE_genes_sCTB, "sCTB", 0);
%sCTB_geneLoc = %$sCTB_geneLoc_ref;
%sCTB_hist = %$sCTB_hist;

#p %aCTB_geneLoc;



#p %DE_genes_sCTB;
# print "$tracklength_aCTB\n";
# print "$tracklength_iCTB\n";
# print "$tracklength_sCTB\n";

my $mirTrackLength = &printMirTrack(0);
&printConnections(0);

&createPhasesFile($tracklength_aCTB, $tracklength_iCTB, $tracklength_sCTB, $mirTrackLength);

system "cat geneNamesTrack.txt mirNamesTrack.txt > namesTrack2.txt";

&printHistogram(\%aCTB_hist, \%iCTB_hist, \%sCTB_hist);

system "tooCircos.sh";

#---------------subs--------------

sub printGeneTracks() {
	# return a data structure that contains the track locations of 
	# genes and mirs.  print gene track and sort by expression
	my $h = shift;
	my %h = %$h;
	my $name = shift;
	my $start = shift;
	my $geneNamesTrack = 'geneNamesTrack.txt';
	my $beast = 'BEAST.txt';
	my %geneLoc = ();  # do this to avoid appending to the hash each time you call the 
						# printGeneTracksFunction. 
	# my $mirNamesTrack = 'mirNamesTrack.txt';
	my %mir_gene_exp = ();
	open(OUT, '>>', $geneNamesTrack);
	open(OUTHIST, '>', $histogramTrack);
	open(BEAST, '>>', $beast);
	# open(OUT3, '>>', $mirNamesTrack);

	# Given the hash made
	# make a hash with mir::target gene as key and it's expression as value
	#print "$name\n";
	#p %h;
	foreach my $i (keys %h) {
		#print "$i\n";
		foreach(@{$h{$i}}) {
		#	print "$_\n";
			(my $gene, my $exp) = split(/\t/, $_);
			$mir_gene_exp{"$i::$gene"} = $exp;
		}
	}
	#p %mir_gene_exp;


	# print "$name\n";
	# p %mir_gene_exp;
	# sort the hash created from above based on exp values. 
	########## I could also try sorting by p-value, this may show a correlation
	########## between number of convergent connections and p-value.  
	#print "sorted\n";
	foreach my $k (sort {$mir_gene_exp{$b} <=> $mir_gene_exp{$a}} keys %mir_gene_exp) { 
		# I could also print histogram files here
		#print "$k\t$mir_gene_exp{$k}\n";	


		my $y = $start + 1;
		(my $mir, my $gene) = split(/::/, $k);
		print BEAST "$name\t$k\t$mir_gene_exp{$k}\n";
		push(@mir6, $mir);
		($alreadyGene{$gene}++) && next;

		print  OUT "$name\t$start\t$y\t0\tname=$gene\n";
		# make a global gene_loc hash
		$geneLoc{"$name::$gene"} = "$start\t$y";
		$geneHist{"$start\t$y"} = $mir_gene_exp{$k};
		$start++;
		
	}
	return ($start, \%geneLoc, \%geneHist);
	close OUT;
	close OUTHIST;
	close BEAST;	
	
	
}
 

sub printMirTrack() {
	# this subroutine takes in all 3 celltype hashes and merges the 
	# mirs together to form the mirtrack
	my $start = shift;
	open(IN, "BEAST.txt");
	my $mirNamesTrack = 'mirNamesTrack.txt';
	open(OUT, '>>', $mirNamesTrack);
	while(<IN>) {
		my ($thing) = ((/\S+/g)[1]);
		# thing actually has all the connections info so I should probably
		# do something with this.  
		my $mir = (split(/::/, $thing))[0];
		push(@mirArray, $mir);
	}
	@mirArray = &removeDups(\@mirArray);
	my $length = @mirArray;
	foreach (sort @mirArray) {

		my $y = $start + 1;
		print OUT "miRNA\t$start\t$y\t0\tname=$_\n";
		# make global mir_loc hash
		$mirLoc{$_} = "$start\t$y";
		$start++;
	}
	close OUT;
	close IN;
	return $length;
}

sub printConnections {
	my $start = shift;
	open(IN, "BEAST.txt");
	my $connections = 'connections.txt';
	open(OUT, ">", $connections);
	system "more BEAST.txt"; # I don't know why, but I need this do go through beast file completely...

	%color = ('aCTB' => 'type=1', 'iCTB' => 'type=2', 'sCTB' => 'type=3');

	while(<IN>) {
		my($thing, $celltype) = ((/\S+/g)[1], (/\S+/g)[0]);
		my($mir, $gene) = (((split(/::/, $thing))[0]), ((split(/::/, $thing))[1]));
		my $celltype_gene = "$celltype::$gene";
		#print "$celltype_gene\n";
			# print "$start\tmiRNA\t$mirLoc{$mir}\t$mir\ttype=1\n";
			# print "$start\t$celltype\t$geneLoc{$celltype_gene}\t$gene\t$color{$celltype}\n";


		###########
		if ($celltype eq "aCTB") { # type 1 = aCTB
			my $celltype_gene = "$celltype::$gene";
			print OUT "$start\tmiRNA\t$mirLoc{$mir}\ttype=1\n";
			print OUT "$start\t$celltype\t$aCTB_geneLoc{$celltype_gene}\ttype=1\n";
		}
		elsif ($celltype eq "iCTB") { # type 2 = iCTB
			my $celltype_gene = "$celltype::$gene";
			print OUT "$start\tmiRNA\t$mirLoc{$mir}\ttype=2\n";
			print OUT "$start\t$celltype\t$iCTB_geneLoc{$celltype_gene}\ttype=2\n";
		}
		elsif ($celltype eq "sCTB") { # celltype eq "sCTB" type 3 = sCTB
			my $celltype_gene = "$celltype::$gene";
			print OUT "$start\tmiRNA\t$mirLoc{$mir}\ttype=3\n";
			print OUT "$start\t$celltype\t$sCTB_geneLoc{$celltype_gene}\ttype=3\n";
		}
		# my $celltype_gene = "$celltype::$gene";
		# print "$celltype\t$mir::$mirLoc{$mir}\t$gene::$geneLoc{$celltype_gene}\n";
		$start++;
	}
	close IN;
	close OUT;
	# p %mirLoc;	
	# p %geneLoc;
}

sub createPhasesFile() {
	my $aCTB = shift;
	my $iCTB = shift;
	my $sCTB = shift;
	my $mirTracklength = shift;
	my $phases = 'phases.txt';
	open(OUTPHASES, ">", $phases);
	print OUTPHASES "chr\t-\taCTB\taCTB\t0\t$aCTB\tgreys-6-seq-3\n";
	print OUTPHASES "chr\t-\tiCTB\tiCTB\t0\t$iCTB\tgreys-6-seq-3\n";
	print OUTPHASES "chr\t-\tsCTB\tsCTB\t0\t$sCTB\tgreys-6-seq-3\n";
	print OUTPHASES "chr\t-\tmiRNA\tmiRNA\t0\t$mirTracklength\tgreys-6-seq-3\n";
	close OUTPHASES;
}

sub printHistogram() {
	my $a = shift;
	my %a = %$a;
	my $i = shift;
	my %i = %$i;
	my $s = shift;
	my %s = %$s;

	my $aCTB_hist = 'aCTB_hist.txt';
	open(OUTA, ">", $aCTB_hist);
	foreach my $x (keys %a) {
		print OUTA "aCTB\t$x\t$a{$x}\n";
	}
	close OUTA;

	my $iCTB_hist = 'iCTB_hist.txt';
	open(OUTI, ">", $iCTB_hist);
	foreach my $y (keys %i) {
		print OUTI "iCTB\t$y\t$i{$y}\n";
	}
	close OUTI;

	my $sCTB_hist = 'sCTB_hist.txt';
	open(OUTS, ">", $sCTB_hist);
	foreach my $z (keys %s) {
		print OUTS "sCTB\t$z\t$s{$z}\n";
	}
	close OUTS;
}

sub removeDups() {
	# fix this for proper dereferencing. 
	my $a = shift;
	my @new = @$a;

	foreach my $p (@new) {
		next if $seen{$p}++;
		push @unique, $p;
	}
	return @unique
}

sub parse_aCTB () {
	my $ctb = shift;
	my %ctb = %$ctb;
	my $mir = shift;
	my %mir = %$mir;
	my %mirGeneExpHash_aCTB;
	foreach my $i (keys %mir) {
		my @compareArray = @{$mir{$i}};
		my $lc = List::Compare->new([keys %ctb], \@compareArray);
		my @intersect = $lc->get_intersection; # get expressed genes with mirs

		if (!@intersect) { # case where mir $i has no target gene DE expressed. 
			next;
		}
		else {
			#print "$i\n";
			# do something with mir $i here.  
			foreach(@intersect) {
				my $gene = $_;
				my $exp = $ctb{$gene};
				my $ge = "$gene\t$exp";
				push(@{$mirGeneExpHash_aCTB{$i}}, $ge);
			}
		}
		# do this later
		# my @Lonly = $lc->get_unique; # get expressed genes w/out mirs
		# my @Ronly = $lc->get_complement; # get NE genes targeted by our mirs.  
		# p @Lonly;	
	}
	#p %mirGeneExpHash;
	# return DE_genes array w/ mirs
	#p @something;
	return \%mirGeneExpHash_aCTB;
}

sub parse_iCTB{
	my $ctb = shift;
	my %ctb = %$ctb;
	my $mir = shift;
	my %mir = %$mir;
	my %mirGeneExpHash_iCTB;
	foreach my $i (keys %mir) {
		my @compareArray = @{$mir{$i}};
		my $lc = List::Compare->new([keys %ctb], \@compareArray);
		my @intersect = $lc->get_intersection; # get expressed genes with mirs

		if (!@intersect) { # case where mir $i has no target gene DE expressed. 
			next;
		}
		else {
			#print "$i\n";
			# do something with mir $i here.  
			foreach(@intersect) {
				my $gene = $_;
				my $exp = $ctb{$gene};
				my $ge = "$gene\t$exp";
				push(@{$mirGeneExpHash_iCTB{$i}}, $ge);
			}
		}
	}
	return \%mirGeneExpHash_iCTB;
}

sub parse_sCTB() {
	my $ctb = shift;
	my %ctb = %$ctb;
	my $mir = shift;
	my %mir = %$mir;
	my %mirGeneExpHash_sCTB;
	foreach my $i (keys %mir) {
		my @compareArray = @{$mir{$i}};
		my $lc = List::Compare->new([keys %ctb], \@compareArray);
		my @intersect = $lc->get_intersection; # get expressed genes with mirs

		if (!@intersect) { # case where mir $i has no target gene DE expressed. 
			next;
		}
		else {
			#print "$i\n";
			# do something with mir $i here.  
			foreach(@intersect) {
				my $gene = $_;
				my $exp = $ctb{$gene};
				my $ge = "$gene\t$exp";
				push(@{$mirGeneExpHash_sCTB{$i}}, $ge);
			}
		}
	}
	return \%mirGeneExpHash_sCTB;
}

# sub createCircosConf() {

# }
######################cool code snippets ###############
# cool sorting:  
# $var = name:2:chrIII
# sort (split(/:/,@var))[N]

########## immediate splitting ############

# my $thing = (split(/\t/, $_))[1];
# my ($thing, $thing2) = ((/\S+/g)[1], (/\S+/g)[0]);
# print "$thing\t$thing2\n"; 
# print

__OUTHIST__
__OUTHIST__
__OUTHIST__
