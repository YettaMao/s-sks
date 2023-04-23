# ---------------------------------------------------
#
# perl script for plotting U V W components 
# 
# Written by Liang Zhao, IGG.CAS
# July 5, 2005
# ---------------------------------------------------

$plot_w = 0;
$offset = 5;

if($plot_w == 1) {
   @comp  = ("u.xy", "w.xy", "v.xy");
   @title = ("x1","x3","x2");
}
else {
   @comp  = ("u.xy","v.xy");
   @title = ("x1","x2");
}

$prem_iaspi = 8;
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
$YMIN = $seg[0]-0.15; 
$YMAX = $seg[1]+0.15;
$XMIN += 5;
$XMAX -= 5;

if($offset>0) {
$XMIN = - $offset;
$XMAX = $XMIN + 50;
}

$sks_t1 = &sks_tab($YMIN);
$sks_t2 = &sks_tab($YMAX);
open(TT,">time.of");
print TT "$sks_t1 $YMIN\n";
print TT "$sks_t2 $YMAX\n";
close(TT);

print "xmin= $XMIN, xmax= $XMAX ; ymin= $YMIN, ymax= $YMAX \n";
$xshift = 3;
$yshift  = 5;
$yanot   = 1;
$xanot   = 50;
if($offset > 0) {$xanot = 20;}

$YLEN       = -12;
$XLEN       = 6.5;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line = "3/80/80/80";
$Time_line = "4/0/0/0t20_20:15";

$textsize = 15;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/20;
$y_offset = ($YMAX - $YMIN)/40;
print "plotting , please wait ...\n";
`gmtset ANNOT_FONT_SIZE_PRIMARY 10`;
for($ii =0; $ii< @comp; $ii++)
{
   if($ii==0) {
         `psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
          &PSTEXT($XMAX + 10, $YMAX + 0.25 ,13, 0, 1, 10, "Time (s)",$PS_file);
          &PSTEXT($XMIN - 11.5, $YMIN +($YMAX - $YMIN)/2,13, 90, 1, 10, "Distance (\260)",$PS_file);
          &PSTEXT($XMIN - 9.5, $YMIN -($YMAX - $YMIN)/20,17, 0, 1, 10, "(b)",$PS_file);
   }
   else 
   {  
       $xshift = $XLEN+0.8; $yshift = 0.;
       $BSCALE = "-B${xanot}/${yanot}wSne";
       `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift -Y$yshift $BSCALE>>$PS_file`;
   }

   print "Now is $title[$ii] ...\n";

   `psxy $comp[$ii] -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
   if($offset <=0) {`psxy time.of -JX$SCALE -R$RSCALE  -N -K -O $PORT -W$Time_line >> $PS_file`;}

   $y = $YMIN + $y_offset;
   &PSTEXT($XMIN+($XMAX-$XMIN)/20,$y ,$textsize, $textangle, 1, 8, $title[$ii],$PS_file);

   `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT $BSCALE>>$PS_file`;
}

`psbasemap -J -R -O $PORT $BSCALE>>$PS_file`;

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

