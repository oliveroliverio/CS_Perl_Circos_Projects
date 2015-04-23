use English;
use warnings;
use Inline::Files;
use Data::Printer;
use List::Compare;
use List::MoreUtils;

system "rm geneNamesTrack.txt";
system "rm pathNamesTrack.txt";

$table = "LMD_table.txt";
open (TABLE1, $table);
# below skips the first line of the file (header)
$firstLine = <TABLE1>;

while(<TABLE1>) {
	$line = $_;
	chomp($line);
	my($pathwayInfo, $genesRef) = &sliceArray(\$line);
	@genesArray = @$genesRef;
	foreach my $z (@genesArray) {
		push(@{$mainHash{$pathwayInfo}}, $z);	
	}	
}
close ($firstLine);

#----------------------------------------------------------------
#		Now I have the mainHash with pathwayInfo as keys, and values
#		as the gene arrays

#			1. Make the gene names track
#---------------------------------------------------------------

#&makeGenesTrack(\%mainHash);
$beastRef = &makeBeast(\%mainHash);
my @beast = @$beastRef;

my ($iCTB_size, $cCTB_size, $sCTB_size, $geneLocRef) = &makeGenesTrack(\@beast);
my %geneLoc = %$geneLocRef;

my ($pathwaySize, $pathwayLocRef) = &makePathTrack();
my %pathwayLoc = %$pathwayLocRef;

&createPhasesFile($iCTB_size, $cCTB_size, $sCTB_size, $pathwaySize);
&makeConnections(\@beast, \%geneLoc, \%pathwayLoc);

system "tooCircos.sh";

#-----------------------------------------------------------
#			----------SUBROUTINES------
#-----------------------------------------------------------
# this sub slices the array to get only the genes to the right of the file. 
sub sliceArray() {
	my $lineRef = shift;
	#dereference $line
	my $line = $$lineRef;
	my @array = ();
	@array = split(/\t/, $line);
	my $pathwayInfo = $array[0] . "::" . $array[1] . "::" . $array[2] . "::" . $array[3];
	# slice the array from the first gene to the very 
	# last item in the array
	@slice = @array[5..$#array];
	#clear slice2 in next line: This was a fucking bitch for some reason
	@slice2 = (); 
	for my $i (0..$#slice) {
		if ($slice[$i] =~ /\S+/) {
			push(@slice2, $slice[$i]);
		}
	}
	return ($pathwayInfo, \@slice2);
	# my($pathway, $celltype, $ratio, $zScore) = ($array[0], $array[1], $array[2], $array[3]);
}

#-------------------------------------------------
#		Make beast hash where keys are the
#		$pathway::$celltype:$ratio and values are
#		genes.  
#			* Input: %mainHash
#
#----------------------------------------------
sub makeBeast() {
	my $h = shift;
	my %h = %$h;
	# accessing values of Hash of arrays
	foreach my $i (keys %h) {
		foreach(@{$h{$i}}) {
			my $gene = $_;
			push(@beastArray, "$i::$gene")
		}
	}
	return (\@beastArray);
}


#---------------------------------------------------
#		Create gene names track for 3 cell types
#			* Input:  %mainHash
#			* Output: textFile and %geneLoc_hash
#---------------------------------------------------
# sub makeGenesTrack_OLD3() {
# 	my $beastRef = shift;
# 	my @beast = @$beastRef;
# 	foreach(@beast) {
# 		my($pathway, $celltype, $blah, $blah2, $gene) = split(/::/, $_);
# 		if ($celltype eq "iCTB") {
# 			push(@geneNameArray_iCTB, $gene);
# 		}
# 		elsif ($celltype eq "cCTB") {
# 			push(@geneNameArray_cCTB, $gene);			
# 		}
# 		elsif ($celltype eq "sCTB") {
# 			push(@geneNameArray_sCTB, $gene);
# 		}
# 	}
# 	@deduped_iCTB = &removeDups(\@geneNameArray_iCTB);
# 	$iCTB_size = @deduped_iCTB;
# 	@sorted_iCTB = sort { $a cmp $b } @deduped_iCTB;

# 	@deduped_cCTB = &removeDups(\@geneNameArray_cCTB);
# 	$cCTB_size = @deduped_cCTB;
# 	@sorted_cCTB = sort { $a cmp $b } @deduped_cCTB;

# 	@deduped_sCTB = &removeDups(\@geneNameArray_sCTB);
# 	$sCTB_size = @deduped_sCTB;
# 	@sorted_sCTB = sort { $a cmp $b } @deduped_sCTB;

# 	p @sorted_sCTB;
# }


sub makeGenesTrack() {
	my $aRef = shift;
	my @a = @$aRef;
	my $geneNamesTrack = 'geneNamesTrack.txt';
	foreach(@a) {
		my @b = split(/::/, $_);
		# b0 pathway, b1 celltype, b4 gene
		if ($b[1] eq "iCTB") {
			push(@iCTB_trackArray, $b[4]);
		}
		elsif ($b[1] eq "cCTB") {
			push(@cCTB_trackArray, $b[4]);
		}
		elsif ($b[1] eq "sCTB") {
			push(@sCTB_trackArray, $b[4]);
		}
	}

	$i_geneCount{$_}++ foreach (@iCTB_trackArray);
	foreach(sort {$i_geneCount{$a} <=> $i_geneCount{$b}} keys %i_geneCount) {
		push(@sorted_iCTB, $_);
	}
	$iCTB_size = @sorted_iCTB;

	$c_geneCount{$_}++ foreach (@cCTB_trackArray);
	foreach(sort {$c_geneCount{$a} <=> $c_geneCount{$b}} keys %c_geneCount) {
		push(@sorted_cCTB, $_);
	}
	$cCTB_size = @sorted_cCTB;
	@sorted_cCTB = reverse (@sorted_cCTB);

	$s_geneCount{$_}++ foreach (@sCTB_trackArray);
	foreach(sort {$s_geneCount{$a} <=> $s_geneCount{$b}} keys %s_geneCount) {
		push(@sorted_sCTB, $_);
	} 
	$sCTB_size = @sorted_sCTB;
	
	

	
	# @deduped_iCTB = &removeDups(iCTB_trackArray);
	# $iCTB_size = @deduped_iCTB;
	# @sorted_iCTB = sort { $a cmp $b } @deduped_iCTB;

	# @deduped_cCTB = &removeDups(cCTB_trackArray);
	# $cCTB_size = @deduped_cCTB;
	# @sorted_cCTB = sort { $a cmp $b } @deduped_cCTB;

	# @deduped_sCTB = &removeDups(sCTB_trackArray);
	# $sCTB_size = @deduped_sCTB;
	# @sorted_sCTB = sort { $a cmp $b } @deduped_sCTB;

	open(OUT, '>>', $geneNamesTrack);
	$start_iCTB = 0;
	foreach(@sorted_iCTB) {
		my $y = $start_iCTB + 1;
		print OUT "iCTB\t$start_iCTB\t$y\t0\tname=$_\n";
		$geneLoc{"iCTB::$_"} = "$start_iCTB\t$y";
		$start_iCTB++;
	}
	$start_cCTB = 0;
	foreach(@sorted_cCTB) {
		my $y = $start_cCTB + 1;
		print OUT "cCTB\t$start_cCTB\t$y\t0\tname=$_\n";
		$geneLoc{"cCTB::$_"} = "$start_cCTB\t$y";
		$start_cCTB++;
	}
	$start_sCTB = 0;
	foreach(@sorted_sCTB) {
		my $y = $start_sCTB + 1;
		print OUT "sCTB\t$start_sCTB\t$y\t0\tname=$_\n";
		$geneLoc{"sCTB::$_"} = "$start_sCTB\t$y";
		$start_sCTB++;
	}
	close OUT;
	return ($iCTB_size, $cCTB_size, $sCTB_size, \%geneLoc);
}



# sub makeGenesTrack_OLD() {
# 	my $h = shift;
# 	my %h = %$h;
# 	#-----------------------
# 	#		Split the mainHash into 3 hashes
# 	#		by cell type
# 	#------------------------
# 	foreach my $i (keys %h) {
# 		my @a = split(/::/, $i);

# 		my $celltype = $a[1];
# 		if ($celltype eq 'iCTB') {
# 			$iCTB_hash{$i} = $h{$i};
# 		}
# 		elsif ($celltype eq 'cCTB') {
# 			$cCTB_hash{$i} = $h{$i};
# 		}
# 		elsif ($celltype eq 'sCTB') {
# 			$sCTB_hash{$i} = $h{$i};
# 		}
# 	}
# 	#----------------------
# 	#		Push all genes into their respective 
# 	#		CTB_arrays
# 	#			* Normally I would then sort
# 	#			* but skip this for now
# 	#---------------------
# 	#-----------------iCTB
# 	foreach my $i (keys %iCTB_hash) {
# 		foreach my $x (@{$iCTB_hash{$i}}) {
# 			push(@iCTB_array, $x);
# 		}
# 	}
# 	@iCTB_array = &removeDups(\@iCTB_array);
# 	my @sorted_iCTB = sort { $a cmp $b } @iCTB_array;

# 	#-----------------sCTB
# 	foreach my $i (keys %sCTB_hash) {
# 		foreach my $x (@{$sCTB_hash{$i}}) {
# 			push(@sCTB_array, $x);
# 		}
# 	}
# 	@sCTB_array = &removeDups(\@sCTB_array);
# 	my @sorted_sCTB = sort { $a cmp $b } @sCTB_array;

# 	#-----------------cCTB
# 	foreach my $i (keys %cCTB_hash) {
# 		foreach my $x (@{$cCTB_hash{$i}}) {
# 			push(@cCTB_array, $x);
# 		}
# 	}
# 	@cCTB_array = &removeDups(\@cCTB_array);
# 	my @sorted_cCTB = sort { $a cmp $b } @cCTB_array;
	

# 	#----------------------------------
# 	#			print to geneNamesTrack.txt
# 	#---------------------------------

# 	open(OUT, '>>', $geneNamesTrack);
# 	open(BEAST, '>>', $beast);
# }

#---------------------------------------------------
#		Create Pathways track sorted by Z-score
#			* input: beasthash
#---------------------------------------------------
sub makePathTrack() {
	#my $aRef = shift;
	#my @a = @$aRef;
	my $pathNamesTrack = 'pathNamesTrack.txt';
	$table2 = "LMD_table.txt";
	open (TABLE2, $table2);
	$firstLine2 = <TABLE2>;
	while(<TABLE2>) {
		my $line = $_;
		chomp($line);
		my @a = split(/\t/, $line);
		(my $pathway) = ($a[0]);
		if ($pathway =~ /\S+/) {
			push(@pathwayArrayPathTrack, $pathway);
		}
	}
	close($firstLine2);
	# --------take care of picking the largest one later-------
	# foreach my $i (sort {$pathwayHash{$a}<=>$pathwayHash{$b}} keys %pathwayHash){
	# 	print "$i\t$pathwayHash{$i}\n";
	# }
	@dedupedPathway = &removeDups(\@pathwayArrayPathTrack);
	@dedupedPathway = reverse @dedupedPathway;
	open (OUT2, ">>", $pathNamesTrack);
	my $start = 0;
	my $size = @dedupedPathway;
	$size = $size+1;
	#foreach my $i (sort {$pathwayHash{$a}<=>$pathwayHash{$b}} keys %pathwayHash){
	foreach(@dedupedPathway){
		my $y = $start+1;
		
		$pathwayLoc{$_} = "$start\t$y";
		print OUT2 "Pathway\t$start\t$y\t0\tname=$_\n";
		$start++;
		# $pathwayHist{"$start\t$y"} = $pathwayHash{$i};
		# make pathway histogram later!!
	}
	close OUT2;
	return($size, \%pathwayLoc);
}

#----------------------------------------------------
#		Create Phases file. 
#----------------------------------------------------
sub createPhasesFile() {
	my $iCTBsize = shift;
	my $cCTBsize = shift;
	my $sCTBsize = shift;
	my $pathwaySize = shift;
	my $phases = 'phases.txt';
	open(OUTPHASES, ">", $phases);
	print OUTPHASES "chr\t-\tiCTB\tiCTB\t0\t$iCTBsize\tgreys-6-seq-3\n";
	print OUTPHASES "chr\t-\tcCTB\tcCTB\t0\t$cCTBsize\tgreys-6-seq-3\n";
	print OUTPHASES "chr\t-\tsCTB\tsCTB\t0\t$sCTBsize\tgreys-6-seq-3\n";
	print OUTPHASES "chr\t-\tPathway\tPathway\t0\t$pathwaySize\tgreys-6-seq-3\n";
	close OUTPHASES;
}

#---------------------------------------------------
#		Make connections file connecting genes to pathways
#		Color links by cell type
#---------------------------------------------------
sub makeConnections() {
	my $beastRef = shift;
	my @beast = @$beastRef;
	my $geneLoc = shift;
	my %geneLoc = %$geneLoc;
	my $pathwayLoc = shift;
	my %pathwayLoc = %$pathwayLoc;
	my $connections = 'connections.txt';

	%color = ('iCTB' => 'type=1', 'cCTB' => 'type=2', 'sCTB' => 'type=3');
	foreach my $i (keys %geneLoc) {
		my ($celltype, $gene) = split(/::/, $i);
		if ($celltype eq "iCTB") {
			$iCTB_con_hash{$gene} = $geneLoc{$i};
		}
		elsif ($celltype eq "cCTB") {
			$cCTB_con_hash{$gene} = $geneLoc{$i};
		}
		elsif ($celltype eq "sCTB") {
			$sCTB_con_hash{$gene} = $geneLoc{$i};			
		}
	}

	open(OUT3, ">", $connections);
	my $conStart = 0;
	foreach(@beast) {
		my @a = split(/::/, $_);
		my ($pathway, $celltype, $gene) = ($a[0], $a[1], $a[4]);
		if ($celltype eq "iCTB") {
			print OUT3 "$conStart\tiCTB\t$iCTB_con_hash{$gene}\ttype=1\n";
			print OUT3 "$conStart\tPathway\t$pathwayLoc{$pathway}\ttype=1\n";
		}
		elsif ($celltype eq "cCTB") {
			print OUT3 "$conStart\tcCTB\t$cCTB_con_hash{$gene}\ttype=2\n";
			print OUT3 "$conStart\tPathway\t$pathwayLoc{$pathway}\ttype=2\n";
		}
		elsif ($celltype eq "sCTB") {
			print OUT3 "$conStart\tsCTB\t$sCTB_con_hash{$gene}\ttype=3\n";
			print OUT3 "$conStart\tPathway\t$pathwayLoc{$pathway}\ttype=3\n";
		}
		$conStart++;
	}
	close OUT3;
}


#--------------------------------
#		2ndary subs
#-----------------------------

sub removeDups() { # and sort
	# fix this for proper dereferencing. 
	@new = ();
	my $a = shift;
	my @new = @$a;

	@unique = (); # I need this to avoid refilling the array. 
	foreach my $p (@new) {
		next if $seen{$p}++;
		push @unique, $p;
	}
	return @unique;
}


## 2do later:
#	* sort genes based on the number of connections they form with pathways
#	* do the same for 
#	* make connection thickness a function of no. of connections. The 
#		more connections, the thicker the link. 
# neurotensin antibody.  
# leptin
