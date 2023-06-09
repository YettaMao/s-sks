@ARGV == 3 || die "Usage plot_seg ps_file $comp1 $comp2\n";
my($ps, $comp1, $comp2) = @ARGV;

$PS_file = $ps;
$phase = "SKS";
$t2 = $comp2;
$comp1 = "*.".$comp1;
$comp2 = "*.".$comp2;

@comp = ($comp1, $comp2);

$read_sac = "/home/zhl/work/s-sks/plot/read_sac";
$sacxy = "sac.xy";
$temp = "tmp.xy";
$pre = 10;
$after = 40;
$norm = 0; 

$XMIN = -$pre;
$XMAX = $after;

$YMIN = 999;
$YMAX = 0;
$bazmin = 180;
$bazmax = 0.;
$max_amp = 0;

   @sacs = `ls $comp[0]`;   
   for($ii=0; $ii<@sacs;$ii++){
       $sac= $sacs[$ii]; chomp($sac);
       @seg = split(/\./, $sac);
       $sta = $seg[0];
       $sac2 = $sta.".".$t2;
       $gcarc = `saclst gcarc <$sac`; chomp($gcarc);
       $baz = `saclst baz <$sac`; chomp($baz);
       $max_v = 0; $max_v2 = 0;
       $depmin = `saclst depmin <$sac`; chomp($depmin); 
       $depmax = `saclst depmax <$sac`; chomp($depmax);        
       if($max_v < &fabs($depmin)) {$max_v = &fabs($depmin);}
       if($max_v < &fabs($depmax)) {$max_v = &fabs($depmax);}
       $depmin = `saclst depmin <$sac2`; chomp($depmin); 
       $depmax = `saclst depmax <$sac2`; chomp($depmax);        
       if($max_v2 < &fabs($depmin)) {$max_v2 = &fabs($depmin);}
       if($max_v2 < &fabs($depmax)) {$max_v2 = &fabs($depmax);}

       if($max_v > $max_v2) { $max_amp[$ii] = $max_v;}
       else {$max_amp[$ii] = $max_v2;}
       $t_r_ration[$ii] = $max_v2/$max_v;

       if($bazmin > $baz) {$bazmin = $baz;}
       if($bazmax < $baz) {$bazmax = $baz;}
       if($ii==0 ) {
           $evdp = `saclst evdp <$sac`; chomp($evdp);
           $nzyear = `saclst nzyear <$sac`; chomp($nzyear);
           $nzjday = `saclst nzjday <$sac`; chomp($nzjday);
           $nzhour = `saclst nzhour <$sac`; chomp($nzhour);
           print "evdp= $evdp\n";
       }
       if($YMIN > $gcarc) {$YMIN = $gcarc;}
       if($YMAX < $gcarc) {$YMAX = $gcarc;}
   }

$title = "Event: ".$nzyear.".".$nzjday.",  baz:"." $bazmin\260 - $bazmax\260";

$YMIN0 = $YMIN;
$YMAX0 = $YMAX;
$YMIN -=0.15;
$YMAX +=0.15;
print "xmin= $XMIN xmax= $XMAX ymin=$YMIN0 ymax= $YMAX0\n";

$xshift = 3;
$yshift = 12;
$yanot  = 1;
$xanot  = 20;

$YLEN   = -10;
$XLEN   = 6;

$PORT   = "-P"; 
$SCALE  = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}Wsne";
$BSCALE2 = "-B${xanot}/${yanot}wsne";
$Line = "4/50/50/50";
$unit = ($YMAX0 - $YMIN0)*2;

for($ii = 0; $ii< @comp; $ii++){
   if($ii==0) {
        `psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE > $PS_file`;
        &PSTEXT($XMAX+($XMAX-$XMIN)/10 ,$YMIN - ($YMAX-$YMIN)/15 ,11, 0, 1, 10, $title,$PS_file);

   }
   else {  
        $xshift = $XLEN +0.3;
        `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift -Y0 $BSCALE2 >>$PS_file`;
#        &PSTEXT($XMIN+($XMAX-$XMIN)/15 ,$YMIN + ($YMAX-$YMIN)/55 ,8, 0, 1, 10, $comp[$ii],$PS_file);
   }
   $num = @sacs;
   for($tt=0; $tt< $num; $tt++) {pop @sacs;}
   @sacs = `ls $comp[$ii]`;
   for($ss=0; $ss<@sacs; $ss++){
      $sac = $sacs[$ss]; chomp($sac);
      @sta = split(/\./,$sac);
      $dt = `saclst delta <$sac`; chomp($dt);
      $gcarc = `saclst gcarc <$sac`; chomp($gcarc);
      $sks_t = &sks_tab($gcarc);
      print "sac: $sac max_v=$max_amp[$ss] sks-t=$sks_t\n";
      $tvsr = substr($t_r_ration[$ss],0,4);
      `$read_sac $sac`;
      &creat_xy_seg($sacxy,$temp,$gcarc, $dt, $sks_t, $unit, $max_amp[$ss],$norm);
#      print "sac: $sac gcarc=$gcarc\n";
      `psxy $temp -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W4/0/128/255t7_7:9  >> $PS_file`;
      if($ii ==1) {
#          &PSTEXT($XMAX+($XMAX-$XMIN)/10 ,$gcarc ,8, 0, 1, 10, "$sta[0]",$PS_file);
#          &PSTEXT($XMAX-($XMAX-$XMIN)/15 ,$gcarc-0.01 ,6, 0, 1, 10, "$tvsr",$PS_file);
      }
      `rm $temp`;
   }
}

sub creat_xy{
   my($sacxy,$temp,$gcarc, $dt, $pre, $after, $unit,$max_amp,$norm) = @_;
   open(INPUT, "$sacxy");
   $line = <INPUT>;
   $num = 0;
   $its_max = 0;
   while($line = <INPUT>){
      chomp($line); 
      @seg = split(" ", $line);
      $time = -$pre + $num * $dt; 
      $amp = $seg[1];
      $n_t[$num] = $time ;
      $n_amp[$num] = $amp;
      if($its_max < &fabs($amp)) {$its_max = &fabs($amp);}
      $num ++;
   }
   close(INPUT);
   open(XY, ">$temp");
   print XY ">>\n";
   for($tt=0; $tt < $num; $tt++){
      if($norm ==0) {$amp = $n_amp[$tt]*$unit/$max_amp + $gcarc;}
      else {$amp = $n_amp[$tt]*$unit/$its_max + $gcarc;}
      printf XY ("%-8.3f %-12.6e\n",$n_t[$tt], $amp);
   }
   close(XY);
}

sub creat_xy_seg{
   my($sacxy,$temp,$gcarc, $dt, $sks_t, $unit,$max_amp,$norm) = @_;
   open(INPUT, "$sacxy");
   $line = <INPUT>;
   $num = 0;
   $its_max = 0;
   while($line = <INPUT>){
      chomp($line); 
      @seg = split(" ", $line);
      $time = -$sks_t + $num * $dt; 
      $amp = $seg[1];
      $n_t[$num] = $time ;
      $n_amp[$num] = $amp;
      if($its_max < &fabs($amp)) {$its_max = &fabs($amp);}
      $num ++;
   }
   close(INPUT);
   open(XY, ">$temp");
   print XY ">>\n";
   for($tt=0; $tt < $num; $tt++){
      if($norm ==0) {$amp = -$n_amp[$tt]*$unit/$max_amp + $gcarc;}
      else {$amp = -$n_amp[$tt]*$unit/$its_max + $gcarc;}
      printf XY ("%-8.3f %-12.6e\n",$n_t[$tt], $amp);
   }
   close(XY);
}

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

sub fabs{
   my($value) = @_;
   if($value > 0) {$abs = $value; }
   else {$abs = - $value;}
   return $abs;
}

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}
