# ---------------------------------------------------
#
# perl script for plotting U V W components 
# 
# Written by Liang Zhao, IGG.CAS
# 20141115
# ---------------------------------------------------

$plot_w = 1;
$offset = 1;

$sks_t_offset = 0;

$column_num = 2; 
if($plot_w == 1) {
   @comp  = ("u.xy", "v.xy", "w.xy");
   @title = ("radial","transverse","vertical");
   $column_num = 3; 
}
else {
   @comp  = ("u.xy","v.xy");
   @title = ("radial","transverse");
}

$prem_iaspi = 0;
if($offset > 0){
   for($cc=0; $cc<@comp;$cc++){
       $ln = 0;
       open(RR,"$comp[$cc]");
       print "$comp[$cc]\n";
       while($lll = <RR>){
            chomp($lll);
            $aline[$ln] = $lll; 
            @segs = split(/ +/,$lll);
            if($segs[1] =~ /^line$/) {
                $n_sac = 1;
                print "$segs[0] $segs[1] $segs[2]\n";
                $gcarc = $segs[2];
            }
            else{
               if($segs[0] >500 && $n_sac == 1) {
                   $sks_t = &sks_tab($gcarc);
                   print "gcarc=$gcarc $segs[0] - $sks_t ,amp=$segs[1]\n";
                   $n_sac = 0;
               }
               if($segs[0] > 500) {$segs[0] = $segs[0] - $sks_t + $prem_iaspi;}
               $aline[$ln]= $segs[0]."  $segs[1]";
            }
            $ln ++;
       }
       close(RR);
       open(OO,">$comp[$cc]");
         for($l=0;$l<$ln;$l++) {print OO "$aline[$l]\n";}
       close(OO);
   }
}

$log   = "log";

$PS_file = "Fig1.ps";

$XMIN    = 0;
$XMAX    = 100;

open(LG, "< log.txt");
$l = <LG>; chomp($l);
@seg  = split(" ",$l);
$XMIN = $seg[0]; $XMAX = $seg[1];
$l = <LG>; chomp($l);
@seg  = split(" ",$l);
$YMIN = $seg[0]-0.2; 
$YMAX = $seg[1]+0.2;
$XMIN += 10;
$XMAX -= 20;

if($offset>0) {
  $XMIN = - 15;
  $XMAX = 45;
}

$sks_t1 = &sks_tab($YMIN)+$sks_t_offset;
$sks_t2 = &sks_tab($YMAX)+$sks_t_offset;
open(TT,">time.of");
print TT "$sks_t1 $YMIN\n";
print TT "$sks_t2 $YMAX\n";
close(TT);

$S_t1 = &S_tab($YMIN)+$sks_t_offset;
$S_t2 = &S_tab($YMAX)+$sks_t_offset;
open(TT,">stime.of");
print TT "$S_t1 $YMIN\n";
print TT "$S_t2 $YMAX\n";
close(TT);

print "xmin= $XMIN, xmax= $XMAX ; ymin= $YMIN, ymax= $YMAX \n";
$xshift = 3;
$yshift  = 5;
$yanot   = 1;
$xanot   = 50;
if($offset > 0) {$xanot = 20;}

$YLEN       = "-15c";
$XLEN       = "6.5c";
if($plot_w == 1) {
   $XLEN       = "4.5c";
}

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line = "3/80/80/250";
$Time_line  = "4/0/0/0t20_20:15";
$STime_line = "4/0/0/255t20_20:15";

$textsize = 15;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/20;
$y_offset = ($YMAX - $YMIN)/40;
print "plotting , please wait ...\n";
`gmtset ANNOT_FONT_SIZE_PRIMARY 10`;
`gmtset MEASURE_UNIT cm`;
for($ii =0; $ii< @comp; $ii++)
{
   if($ii==0) {
         `psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
          &PSTEXT($XMAX + 10, $YMAX + $y_offset*4 ,13, 0, 1, 10, "Time (s)",$PS_file);
          &PSTEXT($XMIN - 4* $x_offset, $YMIN +($YMAX - $YMIN)/2,13, 90, 1, 10, "Distance (\260)",$PS_file);
          &PSTEXT($XMIN - 3.5* $x_offset, $YMIN -($YMAX - $YMIN)/20,17, 0, 1, 10, "(b)",$PS_file);
   }
   else 
   {  
       $xshift = "7.3c"; $yshift = 0.;
       if($plot_w == 1) {
            $xshift = "5.3c"
       }

       $BSCALE = "-B${xanot}/${yanot}wSne";
       `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift -Y$yshift $BSCALE>>$PS_file`;
   }

   print "Now is $title[$ii] ...\n";

   `psxy $comp[$ii] -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
   if($offset <=0) {
       `psxy time.of -JX$SCALE -R$RSCALE  -N -K -O $PORT -W$Time_line >> $PS_file`;
       `psxy stime.of -JX$SCALE -R$RSCALE  -N -K -O $PORT -W$STime_line >> $PS_file`;
   }

   $y = $YMIN + $y_offset;
   &PSTEXT($XMIN+($XMAX-$XMIN)/8,$y ,8, $textangle, 1, 10, $title[$ii],$PS_file);

   `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT $BSCALE>>$PS_file`;
}

`gs $PS_file`;


sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
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

sub S_tab{
   my($gcarc) = @_;
   open(TT,"tcurve.out");
   $l = <TT>; chomp($l);
   @seg = split(" ",$l);
   for($kk=0; $kk<@seg; $kk++){
       if($seg[$kk] =~/^S$/) {$col = $kk;};
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
