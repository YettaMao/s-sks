# perl scripts for drawing surface motion with azimuth
# Written by Liang Zhao, IGG, CAS
# Feb 19, 2007
# Last run on March 1, 2012

$read_sac = "/home/zhl/work/sem/plot/read_sac";

$sacxy = "sac.xy";
$PS_file = "sks.ps";
$draw_model = 1;
$x_o = 8785;
$boundary = 104.8;
$top = 424.86; # based upon semmodel after flatterning
$PFF = ($x_o+462)/111.2; #degree
$sta_num = 46;

# Angle offset
# ang_offset: new - old = 270 - baz -> old = new - 270 + baz
$ang_offset = 0.;

$plot_w = 0;
$offset = 10;

# -----------station names -----------------------
   $res_file = "split_par.txt";
   open(SS,"$res_file");
   $num = 0;
   while($line = <SS>){
      chomp($line);
      @segs = split(/ +/,$line);
      $name[$num]  = substr($segs[0],2,2);
      $gcarcs[$num] = $segs[1];
      $num ++;
   }
   close(SS);

if($plot_w == 1) {
   @comp  = ("u.xy", "w.xy", "v.xy");
   @comp2  = ("u2.xy", "w2.xy", "v2.xy");
   @title = ("x1","x3","x2");
}
else {
   @comp  = ("u.xy","v.xy");
   @comp2  = ("u2.xy","v2.xy");
   @title = ("radial","transverse");
}

$plot_raw = 0;
if($plot_raw==1){
   $pwd = `pwd`; chomp($pwd);
   $plot_seg = "/home/zhl/work/aniso/data/plot_seg2.pl";
   $comp1 = "r"; $comp2 = "t";
#   $rawdir = "/home/zhl/work/aniso/data/velocity/new_v/1998.361.00/";
   $rawdir = "/home/zhl/work/aniso/sws/sws2/cifalps/data3/raw/2013_224_09_49_32/";
   chdir $rawdir;
   print "Now plotting sks phases in $rawdir\n";
   print "perl $plot_seg $PS_file $comp1 $comp2\n";
   `perl $plot_seg $PS_file $comp1 $comp2`;
   `mv $PS_file $pwd`;
   chdir $pwd;
}
else {@comp2  = ("u.xy","v.xy","w.xy");}

$prem_iaspi = 0 ;
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
 #               print "$segs[0] $segs[1] $segs[2]\n";
                $gcarc = $segs[2];
            }
            else{
               if($segs[0] >500 && $n_sac == 1) {
                   $sks_t = &sks_tab($gcarc);
#                   print "gcarc=$gcarc $segs[0] - $sks_t ,amp=$segs[1]\n";
                   $n_sac = 0;
               }
               if($segs[0] > 500) {$segs[0] = $segs[0] - $sks_t + $prem_iaspi;}
               else {$segs[0] = $segs[0] + $prem_iaspi;}
               $aline[$ln]= $segs[0]."  $segs[1]";
            }
            $ln ++;
       }
       close(RR);
       open(OO,">$comp2[$cc]");
         for($l=0;$l<$ln;$l++) {print OO "$aline[$l]\n";}
       close(OO);
   }
}

$log   = "log";

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
$XMIN += 100;
$XMAX -=20;

if($offset>0) {
  $XMIN = - 15;
  $XMAX = 35;
}

$sks_t1 = &sks_tab($YMIN);
$sks_t2 = &sks_tab($YMAX);
open(TT,">time.of");
print TT "$sks_t1 $YMIN\n";
print TT "$sks_t2 $YMAX\n";
close(TT);

print "xmin= $XMIN, xmax= $XMAX ; ymin= $YMIN, ymax= $YMAX \n";

   open(BB,">boundary.xy");
   print BB ">\n";
   print BB "$XMIN $PFF\n";
   print BB "$XMAX $PFF\n";
   close(BB);

$xshift = 3.5;
$yshift  = 12;
$yanot   = 1;
$xanot   = 50;
if($offset > 0) {$xanot = 20;}

$YLEN       = -10;
$XLEN       = 6;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line = "5/40/40/40";
$Time_line = "4/0/0/0t20_20:15";

$textsize = 15;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/20;
$y_offset = ($YMAX - $YMIN)/40;
$x_overlap_offset = -$XLEN - 0.3;
print "plotting , please wait ...\n";
`gmtset ANNOT_FONT_SIZE_PRIMARY 10`;
for($ii =0; $ii< @comp; $ii++)
{
   if($ii==0) {
       if($plot_raw==0) {`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;}
       if($plot_raw==1) {`psbasemap -JX$SCALE -R$RSCALE -K -O $PORT -X$x_overlap_offset $BSCALE>>$PS_file`;}
       &PSTEXT($XMAX + 5, $YMAX + 0.5 ,11, 0, 1, 10, "Time (s)",$PS_file);
       &PSTEXT($XMIN - 11.5, $YMIN +($YMAX - $YMIN)/2,11, 90, 1, 10, "Distance (\260)",$PS_file);
       &PSTEXT($XMIN - 7.5, $YMIN -($YMAX - $YMIN)/20,16, 0, 1, 10, "(b)",$PS_file);
   }
   else 
   {  
       $xshift = $XLEN+0.3; $yshift = 0.;
       $BSCALE = "-B${xanot}/${yanot}wSne";
       `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift -Y$yshift $BSCALE>>$PS_file`;
   }

   print "Now is $title[$ii] ...\n";
   `psxy boundary.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W20/180/180/180 >> $PS_file`;

   `psxy $comp2[$ii] -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
   if($offset <=0) {`psxy time.of -JX$SCALE -R$RSCALE  -N -K -O $PORT -W$Time_line >> $PS_file`;}

   $y = $YMIN + $y_offset;
   &PSTEXT($XMIN+($XMAX-$XMIN)/20,$y ,8, $textangle, 3, 5, $title[$ii],$PS_file);

   if($ii==1){
      for($nj=0; $nj<@name;$nj++){
          $nx = $XMAX + ($XMAX-$XMIN)/10;
          $ny = $gcarcs[$nj];
          &PSTEXT($nx,$ny ,8, 0, 1, 10, $name[$nj],$PS_file);
      }
   }

   `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT $BSCALE>>$PS_file`;
}

if($plot_raw == 1) {`rm @comp2`;}

# ----------------- draw models -------------------------
if($draw_model == 1) {
   print "----- Plotting model -----\n";
   $XMIN = $x_o/111.2+1.5;
   $XMAX = $x_o/111.2 + 1200/111.2 - 2;
   $YMIN = 0;
   $YMAX = 300;

   create_aniso("semmodel",2,"an_zone.xy","an_zone_u.xy",$x_o);

   $xshift = -5.8;
   $yshift = 12.;
   $yanot  = 50;
   $xanot  = 1;
 
   $YLEN   = -11.5/($XMAX-$XMIN)*($YMAX-$YMIN)/111.2/1.2;
   $XLEN   = 11.5;
   $Line = "4/50/50/50";
 
   $PORT   = "-P"; 
   $SCALE  = "$XLEN/$YLEN";
   $RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
   $FSCALE = "-F255/255/255";
   $BSCALE = "-B${xanot}/${yanot}WSne";
   `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift -Y$yshift $BSCALE >>$PS_file`;
   `psxy an_zone.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W8/100/100/100 -G100 >> $PS_file`;
   `psxy an_zone_u.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W5/25/25/25 -G50 >> $PS_file`;
   &PSTEXT($XMIN -1.0 , $YMIN + ($YMAX - $YMIN) /2 ,11, 90, 1, 10, "Depth (km)",$PS_file);
   &PSTEXT($XMIN -0.9 , $YMIN - ($YMAX - $YMIN) /9 ,16, 0, 1, 10, "(a)",$PS_file);
   &PSTEXT($XMIN+1.6 ,$YMIN + 40 ,9, 0, 1, 10, "\370: 135\260 ",$PS_file);
   &PSTEXT($PFF+0.05 ,$YMIN - 50 ,9, 90, 1, 10, "PFF",$PS_file);
   &PSTEXT($XMIN+($XMAX-$XMIN)/2 ,$YMAX + 80 ,11, 0, 1, 10, "Distance (\260)",$PS_file);

   @sacs = `ls *.ce`;
   open(SS, ">sta.xy");
   for($ii=0; $ii<@sacs; $ii++){
       $sac = $sacs[$ii]; chomp($sac);
       $x = `saclst gcarc <$sac`; chomp($x);
       print SS "$x -12\n";
   }
   close(SS);
   `psxy sta.xy -JX$SCALE -R$RSCALE  -N -K -O $PORT -Si0.2 -G0/0/0 >> $PS_file`;
   `psxy arrow.xy -JX$SCALE -R$RSCALE  -N -K -O $PORT -SV0.1c/0.2c/0.13c -G100/100/100 >> $PS_file`;
   `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  $BSCALE >>$PS_file`;
}

# -------------------------------------------------------
# draw split parameters
$plot_sws = 1;
$plot_error = 1;
$raw_split = "/home/zhl/work/aniso/data/sc_result";
#$raw_split = "/home/zhl/work/aniso/data/sc/1998.361.00.txt";
$plot_sc_error = 1;
$res_file = "split_par.txt";

if($plot_sws){
   print "Plotting split parameters\n";
   &plot_split($res_file);
}

sub plot_split{
   my($res_file) = @_;

   $yshift = -18.5;
   $yanot  = 30;
   $xanot  = 1;

   $YMIN = 999; $YMAX = -999;
   $XMIN = 999; $XMAX = -999;

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
   $XMIN -=0.5; $XMAX +=0.2;
   $YMIN = 60;  $YMAX = 180;
   $YLEN   = 4.5;
   $XLEN   = 5.2;
   print "$XMIN $XMAX $YMIN $YMAX\n";

   $PORT       = "-P"; 
   $SCALE      = "$XLEN/$YLEN";
   $RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
   $FSCALE     = "-F255/255/255";
   $BSCALE = "-B${xanot}/${yanot}WSne";
   $Line = "3/30/30/30";

   `psbasemap -JX$SCALE -R$RSCALE -K $PORT -O -Y$yshift $BSCALE>>$PS_file`;
   &PSTEXT($XMIN+($XMAX-$XMIN)/2,$YMIN - ($YMAX - $YMIN)/4.5 ,10, 0, 1, 10, "Distance (\260)",$PS_file);
   &PSTEXT($XMIN-($XMAX-$XMIN)/4.,$YMIN +($YMAX - $YMIN)/2 ,10, 90, 1, 10, "\370 (\260)",$PS_file);

   &PSTEXT($XMIN-($XMAX-$XMIN)/4.2,$YMAX +($YMAX - $YMIN)/20 ,16, 0, 1, 10, "(c)",$PS_file);

# Plot PFF
#   $PFF_line = "10/200/200/200";
#   &PlotBNS($PFF);
#   `psxy temp.xy  -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$PFF_line  >> $PS_file`;

   &read_sws($res_file,1);
   `psxy line.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line  >> $PS_file`;
   `psxy dot.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$Line -Sc0.1  >> $PS_file`;

#   &read_sc($raw_split,1);
#   $fill = "0/128/255";
#   $sc_line = "1/$fill";

#   `psxy line2.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$sc_line  >> $PS_file`;
#   `psxy dot2.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$sc_line -Ss0.11 -G$fill >> $PS_file`;

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

   # Plot PFF
 #  &PlotBNS($PFF);
 #  `psxy temp.xy  -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$PFF_line  >> $PS_file`;

   #&read_sws($res_file,1);
   &read_sws($res_file,0);
   `psxy line.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line  >> $PS_file`;
   `psxy dot.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$Line -Sc0.1  >> $PS_file`;

#   &read_sc($raw_split,0);
#   `psxy line2.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$sc_line  >> $PS_file`;
#   `psxy dot2.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$sc_line -Ss0.11 -G$fill >> $PS_file`;
#   `psbasemap -JX -R -O $PORT  -B >> $PS_file`;
}

sub PlotBNS{
   my($xp) = @_;
   open(BB, ">temp.xy");
   print BB ">\n";
   print BB "$xp $YMAX\n";
   print BB "$xp $YMIN\n";
   close(BB);
}

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
#      print "$ang1 gcarc= $gcarc $phi $dphi\n";
      if($dt>0.5 && $dt< 3.0&& $ddt< 2.0){
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

sub psseis{
   my($offset,$filename,$scale) = @_;
   open(SS, "$filename");
   $line = <SS>; chomp($line);
   @seg = split(" ", $line);
   $num = $seg[1];
   for($kk=0; $kk<$num;$kk++){
      $line = <SS>; chomp($line);
      @seg = split(" ", $line);
      $time[$kk] = $seg[0]; 
      $amp[$kk]  = $offset + $seg[1]*$unit;
   }
   close(SS);
   open(SS, ">$filename");
   print SS "> $num\n";
   for($kk=0; $kk<$num;$kk++){
      print SS "$time[$kk] $amp[$kk]\n"; 
   }
   close(SS);
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

sub create_aniso{
   my($fd_file,$layer,$an_file1,$an_file2,$x_o) = @_;
   #$layer: the anisotropic layer number 

   $back = 1250;
   open(FF, "$fd_file");
   $l = <FF>; chomp($l);
   $l_num = $l;
   for($ll=0;$ll<$l_num;$ll++){ 
      $l  = <FF>; chomp($l);
      $an = <FF>; chomp($l);
      $bn = <FF>; chomp($l);

      if($ll == ($l_num -4)) {
         open(AN,">$an_file1");
         print AN ">\n";
         @an_segs= split(/ +/,$an);
         @bn_segs= split(/ +/,$bn);
         $pp = @an_segs;
         for($li=0;$li<@an_segs-1;$li++){
            if($an_segs[$li] <= $back){
                $lx = ($an_segs[$li]+$x_o)/111.2;
                $ly = $top - $bn_segs[$li];
                print AN "$lx $ly\n";
                if($li==0) {$lx0=$lx; $ly0=$ly;}
            }
            $temp1 = 0;
         }
      }
      if($ll == ($l_num-3)) {
         @an_segs= split(/ +/,$an);
         @bn_segs= split(/ +/,$bn);
         
         for($li=0;$li<@an_segs;$li++){
            $point = @an_segs - 2 - $li;
            $lx = ($an_segs[$point]+$x_o)/111.2;
            $ly = $top - $bn_segs[$point];
            print "";
            if($an_segs[$point] <= $back) { print AN "$lx $ly\n";}
         }
         print AN "$lx0 $ly0\n";
         close(AN);
      }

      if($ll == ($l_num -$layer)) {
         open(AN,">$an_file2");
         print AN ">\n";
         @an_segs= split(/ +/,$an);
         @bn_segs= split(/ +/,$bn);
         $pp = @an_segs;
         for($li=0;$li<@an_segs-1;$li++){
            $lx = ($an_segs[$li]+$x_o)/111.2;
            $ly = $top - $bn_segs[$li];
            print AN "$lx $ly\n";
#            print "$pp $li $lx; $top - $bn_segs[$li] $ly\n";
            if($li==0) {$lx0=$lx; $ly0=$ly;}
         }
      }
      if($ll == ($l_num-1)) {
         @an_segs= split(/ +/,$an);
         @bn_segs= split(/ +/,$bn);
         
         for($li=0;$li<@an_segs-1;$li++){
            $point = @an_segs - 2 - $li;
            $lx = ($an_segs[$point]+$x_o)/111.2;
            $ly = $top - $bn_segs[$point];
            print AN "$lx $ly\n";
#            print "$li $lx $ly\n";
         }
      }
   } 
   print AN "$lx0 $ly0\n";
   close(FF);
   close(AN);
}


