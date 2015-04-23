use English;
use warnings;
use Inline::Files;
use Data::Printer;
use List::Compare;
use List::MoreUtils;

$db = "Homo_sapiens.GRCh38.78.gtf";
%chrHash = ('1' => '1', '2' => '2', '3' => '3', '4' => '4', '5' => '5', '6' => '6', '7' => '7', '8' => '8', '9' => '9', '10' => '10', '11' => '11', '12' => '12', '13' => '13', '14' => '14', '15' => '15', '16' => '16', '17' => '17', '18' => '18', '19' => '19', '20' => '20', '21' => '21', '22' => '22', 'X' => 'X', 'Y' => 'Y');

open(DB, $db);
%c = ();
while(<DB>) {
	$line = $_;
	next if ($line =~ /#/);
	my @a = split(/\t/, $line);
	my($chr, $ensembl, $type, $start, $stop, $stuff) = ($a[0], $a[1], $a[2], $a[3], $a[4], $a[8]);
	if (($type =~ /gene/) && ($ensembl =~ /ensembl/) && (exists $chrHash{$chr})) {
		($gene) = $stuff =~ /gene_name "(.+?)"/;
		#print "$chr\t$gene\t$start\t$stop\n";
		## todo later: how to make 22 arrays including x and y arrays w/ elements being 
		## genes?  
		push(@{$chrHashArray{$chr}}, "$gene\t$start\t$stop");
	}
}

foreach my $i (sort {$a <=> $b} keys %chrHashArray) {
	@sortedArray = sort @{$chrHashArray{$i}};
	foreach(@sortedArray){
		print "$i\t$_\n";
	}
}







#------------------------------
# '1' => '1',
# '2' => '2',
# '3' => '3',
# '4' => '4',
# '5' => '5',
# '6' => '6',
# '7' => '7',
# '8' => '8',
# '9' => '9',
# '10' => '10',
# '11' => '11',
# '12' => '12',
# '13' => '13',
# '14' => '14',
# '15' => '15',
# '16' => '16',
# '17' => '17',
# '18' => '18',
# '19' => '19',
# '20' => '20',
# '21' => '21',
# '22' => '22',
# 'X' => 'X',
# 'Y' => 'Y',





