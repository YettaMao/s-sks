# perl script for creating the initial SEM model with Prem model
# written by Liang Zhao, IGG, CAS
# 20131113
# 20140106: adding function of flattern approximation

# Introduction
# from bottom to top

$flat = 0; # should be consistent with run.par
$prem = "prem_492.par";
$initial = "isomodel";
# nrec: the last layer of GRT
# l_x r_x: left and right boundaries of FD region 
# x_step: step in the x direction
$nrec = 23;
$l_x = -50;
$r_x = 1600;
$x_step = 200;

# (1) read prem.par 
print "1. read prem.par\n";
if(! -f $prem) {
  print "Error: $prem not exist or incrrect directory\n";
}
open(PP,"$prem");
$line=<PP>;
$depth = 0;
$mm = 0;
for($ll=$nrec-2;$ll>=0;$ll--)
{
   $line=<PP>; chomp($line);
   @seg = split(/ +/,$line);
   $vp[$ll] = $seg[1];
   $vs[$ll] = $seg[2];
   $rho[$ll]= $seg[3];
   $thickness[$ll] = $seg[4];
   $th[$mm] = $seg[4];

   $depth +=$thickness[$ll];
   print "layer $ll: vp= $seg[1] vs= $seg[2] depth= $depth mm=$mm \n";
   $mm ++;

}
close(PP);

if($flat ==0){
   $depth = 0;
   for($ii=0;$ii<$nrec -1;$ii++){
        $dep = 0;
       for($jj=0;$jj<$ii;$jj++){
           $dep += $th[$jj];
       }
       $dep += 0.5*$th[$ii];
       $kk = $nrec-2-$ii;
       $f0[$kk] = 6371/(6371-$dep);

       $vp[$kk] *= $f0[$kk];
       $vs[$kk] *= $f0[$kk];
       $rho[$kk]*= $f0[$kk];
       $thickness[$kk] *= $f0[$kk];
       $depth += $thickness[$kk];
       print "No $ii from top to bottom  dep=$dep f0[$kk]= $f0[$kk] thickness=$thickness[$kk] depth=$depth\n";
   }
}

# (2) calculate averge velocity and rho
print "2. calculate averge parameters\n";
# lay: increment of layer number for each layer of FD
$lay = 2;
$iz1 = 0;
$iz2 = $lz1 + $lay;

$num = 0;
$depth = 0;
while($iz1 < $nrec-2){
   $top = $iz2-1;
   print "Average layer: $num $iz1 - $top\n";
   $thick = 0; $tp = 0; $ts=0; $mass=0;
   for($i=$iz1;$i<$iz2;$i++){
       print "Current layer: $i vp=$vp[$i] vs=$vs[$i] rho=$rho[$i]\n";
       $thick += $thickness[$i];
       $tp += $thickness[$i]/$vp[$i];
       $ts += $thickness[$i]/$vs[$i];
       $mass +=$thickness[$i]*$rho[$i];
   }
   $depth += $thick;
   $avp[$num] = $thick/$tp;
   $avs[$num] = $thick/$ts;
   $avrho[$num]=$mass/$thick;
   $dep[$num] = $depth;
   printf("Average: vp=%-8.3f  vs= %-8.3f rho= %-8.3f depth= %-8.3f\n",$avp[$num],$avs[$num],$avrho[$num],$depth);
   $iz1 = $iz2;
   $iz2 = $iz1 + $lay;
   if($iz1>$nrec-2) {
      $iz1= $nrec-2;
      $iz2= $nrec-1;
   }
   $num ++;
}

# (3) output to model file 
$nll = $num;
$np = 0;
$x1 = $l_x;
while($x1<$r_x)
{
    $an[$np] = $x1; 
    $x1 = $x1+$x_step;
#    print "x= $x1\n";
    $np ++;
}
$an[$np] = $r_x; $np++;

open(II,">$initial");
printf II ("%-4d\n",$nll);
for($k=0; $k<$nll;$k++)
{
    printf II ("%-4d %-9.4f %-9.4f %-9.4f\n",$np,$avp[$k],$avs[$k],$avrho[$k]);
    for($j=0;$j<$np;$j++) {printf II ("%-11.4f ",$an[$j]);}
    print II "\n";
    for($j=0;$j<$np;$j++) {printf II ("%-11.4f ",$dep[$k]);}
    print II "\n";
}
close(II);
