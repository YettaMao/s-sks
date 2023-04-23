# perl script for calculate average velocity

$Model = "prem.par";
$depth = 410;

open(PREM,"$Model");
$line=<PREM>;
$thick =0;
$time = 0;
while($thick < $depth){
   $line=<PREM>; chomp($line);
   @seg = split(" ",$line);
   $vp = $seg[0];
   $vs = $seg[1];
   $rho = $seg[2];
   $thick += $seg[3];
   $time += $seg[3]/$vs;
   print "Now depth= $thick vs=$vs time= $time\n";
}
close(PREM);

$average = $thick/$time;
print "Average vs= $average\n";
