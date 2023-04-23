# perl script to create (x) file for receivers

$receiver= "receiver.xy";
$start = 103.5;
$end = 106;

$number = 27;

$step = ($end-$start)/($number - 1);

open(RR,">$receiver");
print RR "$number\n";
for($nn=0;$nn<$number;$nn++){
   $x= $start + $nn*$step;
   printf RR ("%-9.3f ",$x);
}
print RR "\n";
close(RR);
