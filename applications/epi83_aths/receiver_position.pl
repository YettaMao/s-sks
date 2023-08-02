# perl script for generate a file storing receiver position.pl

$xyfile = "receiver.txt"; # using in x-z coordinate system
$poisition= "receiver_long_lat";
$x_o= 5500; # unit: km, being consistent with that in run.par
$pi = 3.1415926536;
$deg_2_dis = 6371*$pi/180;

print "1 degree= $deg_2_dis km\n";

$num = 0;
open(XY,$xyfile);
while($line=<XY>){
   chomp($line);
   @segs=split(/ +/,$line);
   $xxs[$num] = $segs[0];
   $longs[$num] = ($segs[0] + $x_o)/$deg_2_dis;
#   $names[$num] = $segs[1];
   $names[$num] = $num + 1;
   $num ++;
}
close(XY);

open(PP,">$poisition");
print PP "$num\n";
for($ii=0;$ii<@longs;$ii++){
   printf PP ("%-8.3f 0.0   %-6s\n",$longs[$ii], $names[$ii]);
   printf    ("%-8.3f %-8.3f 0.0   %-6s\n",$xxs[$ii], $longs[$ii], $names[$ii]);
}
close(PP);
