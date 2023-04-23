# ------------------------------------------------
# perl script for read prem model and output depth
#
# written by Liang Zhao, IGG, CAS
#
# ------------------------------------------------

$File = "prem.par";
$Depth = "depth.par";

open(PRE, "< $File");
$l = <PRE>; chomp($l);
@seg  = split(" ",$l);
$LAY = $seg[0]; $temp = $seg[1];

for($i=0; $i<$LAY; $i++)
{
   $l = <PRE>; chomp($l);
   @seg  = split(" ",$l);
   $vp[$i] = $seg[0];
   $vs[$i] = $seg[1];
   $rho[$i] = $seg[2];
   $thickness[$i] = $seg[3];
   print "the $i th layer: thickness= $thickness[$i] \n";    
}

close(PRE);

open(DEP, "> $Depth");
$dep = 0;
for($i=0; $i<$LAY; $i++)
{
   $dep = $dep + $thickness[$i];
   $num = $i+1;
   print DEP ("$num layer: depth= $dep vp= $vp[$i] vs= $vs[$i]\n"); 
}
close(DEP);
