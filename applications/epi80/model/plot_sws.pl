# perl scripts for plotting sws results

$res_file = "../split_par.txt";
#$raw_split = "/home/zhl/work/aniso/data/sc_result";
$raw_split = "/home/zhl/work/aniso/data/sc/1998.361.00.txt";
$PS_file = "split.ps";
$plot_error = 0;
$plot_sc_error = 1;
$BNS = 104.8; #degree

# Angle offset
# ang_offset: new - old = 270 - baz -> old = new - 270 + baz
$ang_offset = 15.;

`gmtset ANNOT_FONT_SIZE_PRIMARY 10`;
# 1. plotting fast polarization directions vs gcarc
$xshift = 3;
$yshift = 20;
$yanot  = 30;
$xanot  = 1;

$XMIN = 999; $XMAX = 0;

open(SS,"$res_file");
$num = 0;
while($line = <SS>){
   chomp($line);
   @segs = split(/ +/,$line);
   $name[$num]  = $segs[0];
   $gcarc[$num] = $segs[1];
   $gcarcs[$num]= $segs[1];
   $phi[$num]   = $segs[2];
   $dphi[$num]  = $segs[3];
   $dt[$num]    = $segs[4];
   $ddt[$num]   = $segs[5];
   if($XMIN > $gcarc[$num]) {$XMIN = $gcarc[$num]}
   if($XMAX < $gcarc[$num]) {$XMAX = $gcarc[$num]}
#   print "gcarc = $gcarc[$num]\n";
   $num ++;
}
close(SS);

# ------------- phi: fast polarization -------------------------------------
$XMIN -=0.0; $XMAX +=0.2;
$YMIN = 0;  $YMAX = 120;
$YLEN   = 5;
$XLEN   = 6;
print "$XMIN $XMAX $YMIN $YMAX\n";

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line = "3/30/30/30";

`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
&PSTEXT($XMIN+($XMAX-$XMIN)/2,$YMIN - ($YMAX - $YMIN)/4.5 ,10, 0, 1, 10, "Distance (\260)",$PS_file);
&PSTEXT($XMIN-($XMAX-$XMIN)/5.,$YMIN +($YMAX - $YMIN)/2 ,10, 90, 1, 10, "phi (\260)",$PS_file);

&PSTEXT($XMIN-($XMAX-$XMIN)/4.5,$YMAX +($YMAX - $YMIN)/20 ,13, 0, 1, 10, "(a)",$PS_file);

# Plot BNS
$BNS_line = "10/200/200/200";
&PlotBNS($BNS);
`psxy temp.xy  -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$BNS_line  >> $PS_file`;

&read_sws($res_file,1);
`psxy line.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line  >> $PS_file`;
`psxy dot.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$Line -Sc0.1  >> $PS_file`;

&read_sc($raw_split,1);
$fill = "0/128/255";
$sc_line = "1/$fill";

`psxy line2.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$sc_line  >> $PS_file`;
`psxy dot2.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$sc_line -Ss0.11 -G$fill >> $PS_file`;

# -------------- dt: delay time   -----------------------------------------
$YMIN = 0;  $YMAX = 3.2;

print "$XMIN $XMAX $YMIN $YMAX\n";
$yanot  = 1;
$xanot  = 1;
$xshift = $XLEN*1.3;

$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$BSCALE = "-B${xanot}/${yanot}WSne";
`psbasemap -JX$SCALE -R$RSCALE -K $PORT -O -X$xshift $BSCALE >> $PS_file`;
&PSTEXT($XMIN+($XMAX-$XMIN)/2,$YMIN - ($YMAX - $YMIN)/4.5 ,10, 0, 1, 10, "Distance (\260)",$PS_file);
&PSTEXT($XMIN-($XMAX-$XMIN)/5.5,$YMIN +($YMAX - $YMIN)/2 ,10, 90, 1, 10, "time (s)",$PS_file);

# Plot BNS
&PlotBNS($BNS);
`psxy temp.xy  -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$BNS_line  >> $PS_file`;

#&read_sws($res_file,1);
&read_sws($res_file,0);
`psxy line.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line  >> $PS_file`;
`psxy dot.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$Line -Sc0.1  >> $PS_file`;

&read_sc($raw_split,0);
`psxy line2.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$sc_line  >> $PS_file`;
`psxy dot2.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$sc_line -Ss0.11 -G$fill >> $PS_file`;

sub read_sc{
   my($file,$degree) = @_;
   open(FF,"$file");
   open(LL,">line2.xy");
   open(DD,">dot2.xy");
   while($line = <FF>){
      chomp($line);
      @segs = split(/ +/,$line);
      $sta  = $segs[0];
      $distance = 0;
      $match = 0;
      for($ss=0;$ss<@gcarcs;$ss++){
         if($sta =~/^$name[$ss]$/){
            $distance = $gcarcs[$ss];
            $match = 1;
            print "$sta $name[$ss] $distance \n";
         }
      }
      $phi   = $segs[1];
      if($phi<0) {$phi += 180;}
      $dphi  = $segs[2];
      $dt    = $segs[3];
      $ddt   = $segs[4];
# degree =1: phi; =0: delay time
      if($degree == 0){
         $phi= $dt;
         $dphi = $ddt;
      }
#print "$gcarc $phi $dphi\n";
      if($dt>0.5 && $dt< 3.0 && $match == 1){
         print DD "$distance $phi\n";
         if($plot_sc_error == 1){
            print LL ">\n";
            $yu = $phi + $dphi;
            $yd = $phi - $dphi;
            $xl = $distance - 0.03;
            $xr = $distance + 0.03;
            print LL "$distance $yu\n";
            print LL "$distance $yd\n";
            print LL ">\n";
            print LL "$xl $yu\n";
            print LL "$xr $yu\n";
            print LL ">\n";
            print LL "$xl $yd\n";
            print LL "$xr $yd\n";
         }
      }
   }
   close(DD);
   close(LL);
   close(FF);
}

sub read_sws{
   my($file,$degree) = @_;
   open(FF,"$file");
   open(LL,">line.xy");
   open(DD,">dot.xy");
   while($line = <FF>){
      chomp($line);
      @segs = split(/ +/,$line);
      $name  = $segs[0];
      $gcarc = $segs[1];
      $phi   = $segs[2] + $ang_offset;
      if($phi<0) {$phi += 180;}
      $dphi  = $segs[3];
      $dt    = $segs[4];
      $ddt   = $segs[5];
# degree =1: phi; =0: delay time
      if($degree == 0){
         $phi= $dt;
         $dphi = $ddt;
      }
print "$gcarc $phi $dphi\n";
      if($dt>0.4 && $dt< 3.0&& $ddt< 1.5){
         print DD "$gcarc $phi\n";
         if($plot_error == 1){
            print LL ">\n";
            $yu = $phi + $dphi;
            $yd = $phi - $dphi;
            $xl = $gcarc - 0.03;
            $xr = $gcarc + 0.03;
            print LL "$gcarc $yu\n";
            print LL "$gcarc $yd\n";
            print LL ">\n";
            print LL "$xl $yu\n";
            print LL "$xr $yu\n";
            print LL ">\n";
            print LL "$xl $yd\n";
            print LL "$xr $yd\n";
         }
      }
   }
   close(DD);
   close(LL);
   close(FF);
}

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}

sub PlotBNS{
   my($xp) = @_;
   open(BB, ">temp.xy");
   print BB ">\n";
   print BB "$xp $YMAX\n";
   print BB "$xp $YMIN\n";
   close(BB);
}
