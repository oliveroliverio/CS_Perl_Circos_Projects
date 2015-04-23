use English;
use warnings;
use Inline::Files;
use Data::Printer;
use List::Compare;

$table = "testTable.txt";
open(TABLE1, $table);
while(<TABLE1>){
	$line = $_;
	chomp($line);
	print "$line\n";
}
close TABLE1;

