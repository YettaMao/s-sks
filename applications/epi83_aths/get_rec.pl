# perl script for obtaining receiver's positions

$record = "/home/maoyt/work/s-sks/ncisp6/data/raw/2007_289_21_05_41/";
$rec = "receiver.xy";
$sta = "sta_name";
$xmin=9000; # km must be consistent with that in run.par
$degree_to_km = 6371*3.14159265358979/180.; 
$change_to_km = 1;

open(RR,">$rec");
open(SS,">$sta");

chdir $record;
@bhn = `ls *.r`;
$num = @bhn;

# calculating gcarc
for($ii=0;$ii<@bhn;$ii++){
   $sac = $bhn[$ii]; chomp($sac);
   my($sac,$gcarc) = split /\s+/, `saclst gcarc f $sac`; 
   $rec[$ii] = $gcarc;
   @segs = split(/\./,$sac); 
   $name[$ii] =$segs[0]; 
#   print "$sac 's gcarc= $gcarc\n";
}

#
print "sort... \n";
for($ii=$num-1;$ii>=0;$ii--){
   for($jj=0;$jj<=$ii;$jj++){
      if($rec[jj] > $rec[$ii]) {
          $tt = $rec[$ii];
          $tname = $name[$ii];
          $rec[$ii] = $rec[$jj];
          $rec[$jj] = $tt;
          $name[$ii] = $name[$jj];
          $name[$jj] = $tname;
      }
   }
}

printf RR ("%-4d\n",$num);
for($ii=0;$ii<$num;$ii++){
   $sac = $bhn[$ii]; chomp($sac);
   $xx = $rec[$ii];
   if($change_to_km ==1){
      $xx = $xx * $degree_to_km - $xmin;
   }
   printf RR ("%-7.3f ",$xx);
   print SS "$name[$ii]\n";
   print "No.$ii $sac: $rec[$ii] , $xx\n";
}
print RR "\n";

close(RR);
close(SS);
