# perl script to create (x) file for receivers

$receiver= "receiver.xy";
$start = 82.4;
$end = 90.8;

$number = 60;

$step = ($end-$start)/($number - 1);

open(RR,">$receiver");
print RR "$number\n";
for($nn=0;$nn<$number;$nn++){
   $x= $start + $nn*$step;
   $x = $x*111.2-9000;
   printf RR ("%-9.3f ",$x);
}
print RR "\n";
close(RR);
