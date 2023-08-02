@ARGV == 3 || die "Usage plot_seg ps_file $comp1 $comp2\n";
my($ps, $comp1, $comp2) = @ARGV;
#use strict;
#use warnings;

$PS_file = $ps;
$phase = "SKS";
$t2 = $comp2;
$comp1 = "*.".$comp1;
$comp2 = "*.".$comp2;

@comp = ($comp1, $comp2);

$read_sac = "/home/maoyt/work/s-sks/plot/read_sac";
$read_sac = "/home/maoyt/work/s-sks/plot/readsac_seg";
$sacxy = "sac.xy";
$temp = "tmp.xy";
$pre = 10;
$after = 50;
$norm = 0; 

$XMIN = -$pre;
$XMAX = $after;

$YMIN = 999;
$YMAX = 0;
$bazmin = 180;
$bazmax = 0.;
$max_amp = 0;
$ii = 0;

@sacs = `ls $comp[0]`;   
for($ii=0; $ii<@sacs;$ii++){
    $sac= $sacs[$ii]; chomp($sac);
    @seg = split(/\./, $sac);
    $sta = $seg[0];
    $sac2 = $sta.".".$t2;
    ($sac,$gcarc,$baz) = split /\s+/,`saclst gcarc baz f $sac`; 
    $sks_t = &sks_tab($gcarc);
    $t_b = $sks_t - $pre;
    $t_e = $sks_t + $after;
    `$read_sac $sac $t_b $t_e`; 
    $max_v = &readmax() ;chomp($max_v);
    `$read_sac $sac2 $t_b $t_e`; 
    $max_v2=&readmax() ; chomp($max_v2);
    print "No. $ii max_v=$max_v v2=$max_v2\n";
    if($max_v > $max_v2) { $max_amp[$ii] = $max_v;}
    else {$max_amp[$ii] = $max_v2;}
    $t_r_ration[$ii] = $max_v2/$max_v;

    if($bazmin > $baz) {$bazmin = $baz;}
    if($bazmax < $baz) {$bazmax = $baz;}
    if($ii==0 ) {
        ($sac,$evdp,$nzyear,$nzjday,$nzhour) = split /\s+/,`saclst evdp nzyear nzjday nzhour f $sac`; 
        print "evdp= $evdp\n";
    }
    if($YMIN > $gcarc) {$YMIN = $gcarc;}
    if($YMAX < $gcarc) {$YMAX = $gcarc;}
}

$title = "Event: ".$nzyear.".".$nzjday.",  baz:"." $bazmin\260 - $bazmax\260";

$YMIN0 = $YMIN;
$YMAX0 = $YMAX;
$YMIN -=0.25;
$YMAX +=0.25;
print "xmin= $XMIN xmax= $XMAX ymin=$YMIN0 ymax= $YMAX0\n";

$xshift = 3.5;
$yshift = 4;
$yanot  = 1;
$xanot  = 15;

$YLEN   = -15;
$XLEN   = 6;

$PORT   = "-P"; 
$SCALE  = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE = "-F255/255/255";
#$BSCALE = "-B${xanot}/${yanot}wsne";
#$BSCALE2 = "-B${xanot}/${yanot}wsne";
$BSCALE = "-BWSne";
$BSCALE2 = "-BwSne";
$Line = "4/50/50/50";
$unit = ($YMAX0 - $YMIN0)/20;

`gmt begin sks ps`;
for($ii = 0; $ii< @comp; $ii++){
   if($ii==0) {
        `gmt basemap -JX$SCALE -R$RSCALE -X$xshift -Y$yshift $BSCALE -Bxa15f15 -Bya1`;
        &PSTEXT($XMAX+($XMAX-$XMIN)/10 ,$YMIN - ($YMAX-$YMIN)/15 ,11, 0, 1, 10, $title,$PS_file);

   }
   else {  
        $xshift = $XLEN +0.3;
        `gmt basemap -JX$SCALE -R$RSCALE -X$xshift -Y0 $BSCALE2 -Bxa15f15 `;
	# &PSTEXT($XMIN+($XMAX-$XMIN)/15 ,$YMIN + ($YMAX-$YMIN)/55 ,8, 0, 1, 10, $comp[$ii],$PS_file);
   }
   $num = @sacs;
   for($tt=0; $tt< $num; $tt++) {pop @sacs;}
   @sacs = `ls $comp[$ii]`;
   for($ss=0; $ss<@sacs; $ss++){
      $sac = $sacs[$ss]; chomp($sac);
      @sta = split(/\./,$sac);
      ($sac,$dt,$gcarc) = split /\s+/,`saclst delta gcarc f $sac`; 
      $sks_t = &sks_tab($gcarc);
      $t_b = $sks_t - $pre;
      $t_e = $sks_t + $after;

      # print "sac: $sac max_v=$max_amp[$ss]  sks-t=$sks_t\n";
      $tvsr = substr($t_r_ration[$ss],0,4);
      `$read_sac $sac $t_b $t_e`;
      &creat_xy_seg($sacxy,$temp,$gcarc, $dt, $sks_t, $unit, $max_amp[$ss],$norm);

      print "sac: $sac gcarc=$gcarc\n";
      `gmt plot $temp -JX$SCALE -R$RSCALE -W1.2p,deepskyblue,: -N `;
      if($ii ==1) {
          &PSTEXT($XMAX+($XMAX-$XMIN)/10 ,$gcarc ,8, 0, 1, 10, "$sta[0]",$PS_file);
	  #  &PSTEXT($XMAX-($XMAX-$XMIN)/15 ,$gcarc-0.01 ,6, 0, 1, 10, "$tvsr",$PS_file);
      }
#      `rm $temp`;
   }
}
`gmt end`;
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
      $time = -$sks_t + $seg[0]; 
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
   open(MYT, ">>obs_t.xy");
   print MYT ">>\n";
   for($tt=0; $tt < $num; $tt++){
      if($norm ==0) {$amp = -$n_amp[$tt]*$unit/$max_amp + $gcarc;}
      else {$amp = -$n_amp[$tt]*$unit/$its_max + $gcarc;}
      printf MYT ("%-8.3f %-12.6e\n",$n_t[$tt], $amp);
   }
   close(MYT);
}

sub sks_tab{
   my($gcarc) = @_;
   open(TT,"taup_curve.gmt");
   $l = <TT>; chomp($l);
   @seg = split(" ",$l);
   $num = 0;
   while($l=<TT>){
      chomp($l);
      @seg = split(" ",$l);
      $g[$num] = $seg[0];
      $sks[$num] = $seg[1];
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
    open(GMT,"| gmt text -R$RSCALE -N -JX$SCALE -F+f+a+j ");
        print GMT "$xx $yy $textsize,$textfont $textangle $just $text\n";
    close(GMT);
}

sub readmax{
   open(MAX, "maxv.txt");
   $line=<MAX>;
   chomp($line);
   return $line;
   close(MAX);
}
