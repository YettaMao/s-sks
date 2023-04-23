# perl scropt for transferring isotropic model to anisotropic model
# written by Liang Zhao, IGG, CAS
# 20131113

$iso_model = "isorf";
$aniso_model = "semmodel";

open(II, "$iso_model");
open(AA, ">$aniso_model");
$line = <II>; chomp($line);
$nll = $line;
printf AA ("%-4d\n",$nll);

for($ll=0;$ll<$nll;$ll++){
   $line = <II>; chomp($line);
   @seg = split(/ +/, $line);
   $np = $seg[0];
   $vp = $seg[1]; 
   $vs = $seg[2];
   $rho= $seg[3];
   $c11= $vp*$vp*$rho;
   $c16= 0.;
   $c22= $c11;
   $c26= 0.0;
   $c33= $c11;
   $c36= 0.;
   $c44= $vs*$vs*$rho;
   $c45= 0.0;
   $c12= $c11 - 2*$c44;
   $c13= $c11 - 2*$c44;
   $c23= $c11 - 2*$c44;
   $c55= $c44;
   $c66= $c44;

   printf AA ("%-4d %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f",$np,$c11,$c12,$c13,$c16,$c22,$c23,$c26);
   printf AA ("%-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f %-10.5f\n",$c33,$c36,$c44,$c45,$c55,$c66,$rho);

   $line = <II>; chomp($line);
   @seg = split(/ +/, $line);
   for($kk=0;$kk<@seg;$kk++){printf AA ("%-11.4f ",$seg[$kk]);}
   print AA "\n";

   $line = <II>; chomp($line);
   @seg = split(/ +/, $line);
   for($kk=0;$kk<@seg;$kk++){
       $y = $seg[$kk] - 0.;
       printf AA ("%-11.4f ",$y);
   }
   print AA "\n";
}

close(II);
close(AA);
