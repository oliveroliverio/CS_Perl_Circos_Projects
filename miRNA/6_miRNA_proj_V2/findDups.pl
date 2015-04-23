use English; 
use warnings;

$file = "DE_expression_3_celltypes.txt";
$/ = ">";
open(IN, $file);
while(<IN>) {
	chomp;
	/\n/;
	$header = $PREMATCH;
	$data = $POSTMATCH;
	($celltype, $m) = split(/\t/, $header); 
	print "$celltype\n";
}
