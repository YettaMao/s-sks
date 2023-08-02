# perl for calculate splitting parameters of out put files
# call Silver and Chan's method [1991]
# By Liang Zhao, IGG, CAS 
# Nov 2, 2006

$split_par = "split_par.txt";
$SC = "/home/zhl/work/aniso/sws/local/bin/an_split";
$before = 10;
$after = 30;

$log = "sc_result";
$name_file = "sta_name";

open(NN,"$name_file");
$num = 0;
while($line = <NN>){
   chomp($line);
   $name[$num] = $line;
   $num ++;
}
close(NN);

@sacs = `ls *.ce`;

if(@sacs != $num) {print "Error: station number is not equal to receivers\n";}

open(SS, ">$split_par");
for($ii=0; $ii<@sacs; $ii++){
#   $sac = $sacs[$ii]; chomp($sac);
#   @seg = split(/\./,$sac); 
#   $sta = $seg[0];
   $sta = $ii;
   $sac = $sta.".ce";
   $gcarc = `saclst gcarc <$sac`; chomp($gcarc);
   $omark = `saclst o <$sac`; chomp($omark);
   $sks_t = &sks_tab($gcarc);
   $n_sks = $sks_t + $omark;
   $n_sks = $sks_t;
#   print "$sta gcarc=$gcarc, O= $omark sks arrive at $sks_t -> $n_sks\n";
   $A = $n_sks - $before;
   $F = $n_sks + $after;
   print "call pgs: $SC $sta $A $F\n";
   `$SC $sta $A $F`;

# read splitting parameters
   open(LL, "$log");
   $line =<LL>;
   chomp($line);
   @seg = split(/ +/, $line);
      $phi  = $seg[0];
      $d_phi= $seg[1];
      $dt   = $seg[2];
      $d_dt = $seg[3];
   close(LL);
   print "$sta phi=$phi d_phi=$d_phi dt=$dt d_dt=$d_dt\n";
   printf SS ("%-5s %-7.2f %-6.1f %-6.1f %-7.2f %-7.2f\n", $name[$ii], $gcarc, $phi, $d_phi, $dt, $d_dt);
}
close(SS);

#`rm $log`;

sub sks_tab{
   my($gcarc) = @_;
   open(TT,"tcurve.out");
   $l = <TT>; chomp($l);
   @seg = split(" ",$l);
   for($kk=0; $kk<@seg; $kk++){
       if($seg[$kk] =~/^SKSac$/) {$col = $kk;};
   }
#   print "sks at $col column\n";
   $num = 0;
   while($l=<TT>){
      chomp($l);
      @seg = split(" ",$l);
      $g[$num] = $seg[0];
      $sks[$num] = $seg[$col];
#      print "$g[$num] $sks[$num]\n";
      $num ++;
   }
   for($kk=0; $kk<@g-1;$kk++){ 
      if($gcarc>=$g[$kk] && $gcarc< $g[$kk+1]){
          $time = $sks[$kk] + ($gcarc-$g[$kk])/($g[$kk+1]-$g[$kk])*($sks[$kk+1] - $sks[$kk]);
      }
   }
   close(TT);
   return $time;
}
