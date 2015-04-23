# need to take in 2 files:  connections file and file
# that associates miR name and location.  if miRname is 
# any of these: mir-21, and so on.  then color it.  
# if location matches, then put 
$connections = $ARGV[0];
open(IN1, $connections);
while(<IN1>) {
	chomp;
	$line = $_;
	print "$line\n";
}
