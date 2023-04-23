# perl script for bat run calculation program
# step 1: creating new cijkl matrix by invoke elastic_rotate
#
#
# Last Modify: zhl, Oct 20, 2006

$m_dir = "/home/zhl/work/sem/model/src/";
$rot = "/home/zhl/work/sem/model/src/elastic_rotate";
$model = "/home/zhl/work/sem/model/src/TI.ela";
$new_model = "new.ela";

$sem = "semmodel_old";
$sem_new = "semmodel";

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

$work = ` pwd`; chomp($work);

@angs = (0,5,10,15,20,25,30,35,40,45,50,55,60,65,68,70,75,80,85,175);
for($aa=0; $aa<@angs; $aa++){
# step 1 -------------------------------
   $angle_100 = $angs[$aa];
   chdir $m_dir;
# $angle is initialy refter the rotating angle of the new coordinate system WRT old system
   $angle = $angle_100 - 90;
   $sem_new = "sem_".$angle_100;

   print "\n No.$aa rot wrt old system angle=$angle, wrt North= $angle_100\n";
   print " $rot $model $new_model $angle \n";
   `$rot $model $new_model $angle`;
   open(MM, "$new_model");
   for($i=0;$i<6;$i++){
       $line=<MM>; chomp($line);
       @seg = split(" ",$line);
       if(@seg != 6) {print "wrong elastic matrix\n";}
       for($j=0;$j<@seg;$j++){
          $cij[$i*6+$j] = $seg[$j];
#          print "$cij[$i*6+$j] "; 
       }
#       print "\n";
   }
   close(MM);
   print "\n";

# step 2 -------------------------------
   chdir $work;
   print "compile new SEM model: sem\n";
   &new_sem($sem,$sem_new,@cij);

}

sub new_sem{
   my($file,$newfile,@ckl) = @_;
   open(FF, "$file");
   $l = <FF>; chomp($l);
   $layer = $l;
   for($ll=0;$ll<$layer;$ll++){ 
      $l = <FF>; chomp($l);
      $lines[$ll*3] = $l;
      if($ll == $layer -1) {
         @ss = split(" ", $l);
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
         $newl = $ss[0];
         for($si=1;$si<@ss;$si++){
            $newl = $newl."    ".$ss[$si];
         }
         $lines[$ll*3] = $newl;
         print "\n $newl \n";
      }
      $l = <FF>; chomp($l);
      $lines[$ll*3+1] = $l;
      $l = <FF>; chomp($l);
      $lines[$ll*3+2] = $l;
   } 
   close(FF);
   print "open $newfile\n";
   open(FF,">$newfile");
   print FF "$layer\n";
   for($si=0;$si<@lines;$si++){
      print FF "$lines[$si]\n";
   }
   close(FF);
}

