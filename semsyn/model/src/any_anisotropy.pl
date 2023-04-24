# perl script for generating any anisotropy medium from an isotropic medium
# and calculating velocity along different direction of propagation
# By Liang Zhao, IGG, CAS
# cp cal_v.pl
# Ref: Crampin, S.,1977, A review of the effects of anisotropic layering on the propagation of seismic waves, Geophys. J. R. astr. Soc., 49, 9-27., Page10
# 20140319 

# consider plane wave propagating along the xn axis: uj=aj exp[iw(t-pk xk/c)], j,k = 1,2,3
# eigenvalue proble: (Cinjn - lama*E)a = 0 (i,j=1,2,3)
$cal = "./cjcbi0";
$n=3; 
$stronger = 1.3; # enlarge at 1 direction (x)

#@vps =(7.30);
#@vss =(4.34);
#@rhos=(3.00);

#@vps =(7.80);
#@vss =(4.50);
#@rhos=(3.30);

#@vps =(6.90);
#@vss =(3.96);
#@rhos=(2.0);

@vps =(8.10);
#@vss =(4.73);
@vss =(4.48);
@rhos=(3.2);

$TI = "any.ela";

$rot = "/home/maoyt/work/rf/semsyn/model/src/elastic_rotate";
$new_model = "newany.ela";

$angle_100 = 175; # fast polarization direction wrt north 
$sem = "sem_".$angle_100;

$work = ` pwd`; chomp($work);

###############################
#     ^ North
#     |
#     |
#     |
#     |
#     |
#     |------------------> x or east
#  elastic_rotate: angle is wrt to x
#  then wrt North, angle += 90

for($ii=0;$ii<@vps;$ii++){
   $vs  = $vss[$ii];
   $vp  = $vps[$ii];
   $rho = $rhos[$ii];

#  step 0: generate cij matrix to $model
   $c11= $vp*$vp*$rho*$stronger;
   $c16= 0.;
   $c22= $vp*$vp*$rho;
   $c26= 0.0;
   $c33= $vp*$vp*$rho;
   $c36= 0.;
   $c44= $vs*$vs*$rho*$stronger;
   $c45= 0.0;
   $c12= $c11 - 2*$c44;
   $c13= $c11 - 2*$c44;
   $zero = 0.;

   $c55= $vs*$vs*$rho;
   $c66= $vs*$vs*$rho;
   $c23= $c22 - 2*$c55;

   print "   c11     c12      c13      c16      c22      c23      c26     c33      c36      c44      c45      c55      c66     rho\n";
   printf ("%-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f",$c11,$c12,$c13,$c16,$c22,$c23,$c26);
   printf ("%-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f\n",$c33,$c36,$c44,$c45,$c55,$c66,$rho);


   open(OT,">$TI");
   printf OT ("%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",$c11,$c12,$c13,$zero,$zero,$zero);
   printf OT ("%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",$c12,$c22,$c23,$zero,$zero,$zero); 
   printf OT ("%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",$c13,$c23,$c33,$zero,$zero,$zero);
   printf OT ("%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",$zero,$zero,$zero,$c44,$zero,$zero);
   printf OT ("%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",$zero,$zero,$zero,$zero,$c55,$zero);    
   printf OT ("%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f\n",$zero,$zero,$zero,$zero,$zero,$c66);        
   close(OT);

# step 1: calculate anisotropic parameters
# Usage: ./cjcbi0 TI.ela  3 3.3
   print "$cal $TI $n $rho\n";
`  $cal $TI $n $rho`;

# step 2: rotate the TI matrix

   $angle = $angle_100 - 90;
   print "rot wrt old system angle=$angle, wrt North= $angle_100\n";
   print "$rot $TI $new_model $angle\n";
   `$rot $TI $new_model $angle`;
   open(MM, "$new_model");
   for($i=0;$i<6;$i++){
       $line=<MM>; chomp($line);
       @seg = split(" ",$line);
       if(@seg != 6) {print "wrong elastic matrix\n";}
       for($j=0;$j<@seg;$j++){
          $cij[$i*6+$j] = $seg[$j];
          printf ("%-9.3f ",$cij[$i*6+$j]); 
       }
       print "\n";
   }
   close(MM);
   print "\n";

# step 3: output to model file
   print "compile new SEM model: sem\n";
   &new_sem($vp,$vs,$rho,$sem,@cij);
}


sub new_sem{
   my($vp,$vs,$rho,$newfile,@ckl) = @_;

         $ss[1]  = $ckl[0*6+0];
         $ss[2]  = $ckl[0*6+1];
         $ss[3]  = $ckl[0*6+2];
         $ss[4]  = $ckl[0*6+5];
         $ss[5]  = $ckl[1*6+1];
         $ss[6]  = $ckl[1*6+2];
         $ss[7]  = $ckl[1*6+5];
         $ss[8]  = $ckl[2*6+2];
         $ss[9]  = $ckl[2*6+5];
         $ss[10] = $ckl[3*6+3];
         $ss[11] = $ckl[3*6+4];
         $ss[12] = $ckl[4*6+4];
         $ss[13] = $ckl[5*6+5];
         $newl="";
         for($si=1;$si<@ss;$si++){
            $newl = $newl.$ss[$si]."    ";
         }
         $newl = $newl."    ".$rho;
   print "open $newfile\n";
   open(FF,">$newfile");
   print FF "$vp $vs $rho\n";
   print FF "$newl\n";
   print "$vp $vs $rho\n";
   print "   c11     c12      c13      c16      c22      c23      c26     c33      c36      c44      c45      c55      c66       rho\n";
   print "$newl\n";
   close(FF);
}
