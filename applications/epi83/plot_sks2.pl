# perl scripts for drawing surface motion with azimuth
# Written by Liang Zhao, IGG, CAS
# Feb 19, 2007
# Last run on 20141117

$read_sac = "/home/zhl/work/sem/plot/read_sac";

$sacxy = "sac.xy";
$PS_file = "sks2.ps";
$draw_model = 1;
$x_o = 9000;
$boundary = 104.8;
$top = 446.2710; # based upon semmodel after flatterning
$PFF = ($x_o+462)/111.2; #degree
$sta_num = 52;
$pi = 3.14159265;
$deg = 6371*$pi/180.;

$offset_S = 0.8; # due to the low verlocity above the CMB
$offset_S_b = 0.02; # due to the low verlocity above the CMB

# Angle offset
# ang_offset: new - old = 270 - baz -> old = new - 270 + baz
$ang_offset = 0.;

$plot_w = 0;
$offset = 5;
$prem_iaspi = 3 ;

`cp ../vel_gray.cpt ./`;

# -----------station names -----------------------
   $res_file = "split_par.txt";
   open(SS,"$res_file");
   $num = 0;
   while($line = <SS>){
      chomp($line);
      @segs = split(/ +/,$line);
      $name[$num]  = substr($segs[0],3,2);
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
   @comp3  = ("u3.xy", "w3.xy", "v3.xy");
   @title = ("radial","transverse");
}

$plot_raw = 1;
if($plot_raw==1){
   $pwd = `pwd`; chomp($pwd);
   $plot_seg = "/home/zhl/work/s-sks/plot/plot_seg3.pl";
   $comp1 = "r"; $comp2 = "t";
   $rawdir = "/home/zhl/work/s-sks/ncisp6/data/raw/2007_289_21_05_41/";
   chdir $rawdir;
   print "Now plotting sks phases in $rawdir\n";
   print "perl $plot_seg $PS_file $comp1 $comp2\n";
   `perl $plot_seg $PS_file $comp1 $comp2`;
   `mv $PS_file $pwd`;
   chdir $pwd;
}
else {@comp2  = ("u.xy","v.xy","w.xy");}

if($offset > 0){

   for($cc=0; $cc<@comp;$cc++){
       $time_cut = 8;
       open(OO,">$comp2[$cc]");
       open(OFS,">$comp3[$cc]");
       $ln = 0;
       $ln1=0;
       open(RR,"$comp[$cc]");
       print "$comp[$cc]\n";
       while($lll = <RR>){
            chomp($lll);
            $aline[$ln] = $lll; 
            @segs = split(/ +/,$lll);
            $sks_t =0;
            if($segs[1] =~ /^line$/) {
#                print "$segs[1]\n";
                $n_sac = 1;
                $gcarc = $segs[2];
                print OO ">> \n";
                print OFS ">> \n";
            }
            else{
               if($segs[0] >500 && $n_sac == 1) {
                   $sks_t = &sks_tab($gcarc);
#                   print "gcarc=$gcarc $segs[0] - $sks_t ,amp=$segs[1]\n";
                   $n_sac = 0;
               }
               $segs[0] = $segs[0] - $sks_t + $prem_iaspi;
               if($gcarc < 84) {$time_cut = 6;}
               else {   $time_cut = 8;}
#  after the 8 s, it is S
               if($segs[0] >= $time_cut){
                  $goff= $gcarc - 82;
                  $segs[0] = $segs[0] + $goff * $offset_S + $goff * $goff * $offset_S_b;
                  $aline1= $segs[0]."  $segs[1]";
                  print OFS "$aline1\n";
                  $ln1 ++;
               }
               else{
                  $aline= $segs[0]."  $segs[1]";
                  print OO "$aline\n";
                  $ln ++;
               }
# add offset caused by low velocity above the CMB
           }
       }
       close(RR);

       close(OO);
       close(OFS);
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
$YMIN = $seg[0]; 
$YMAX = $seg[1];
$YMIN -=0.25;
$YMAX +=0.25;

$XMIN += 100;
$XMAX -=20;

if($offset>0) {
  $XMIN = - 10;
  $XMAX = 12;
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
$yshift  = 4;
$yanot   = 1;
$xanot   = 5;
if($offset > 0) {$xanot = 5;}

$YLEN       = -15;
$XLEN       = 6;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line = "3/80/80/80";
$Line = "4/200/30/30";
$Liner = "5/200/40/40";
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
       &PSTEXT($XMIN - 3.5, $YMIN -($YMAX - $YMIN)/20,15, 0, 1, 10, "(b)",$PS_file);
   }
   else 
   {  
       $xshift = $XLEN+0.3; $yshift = 0.;
       $BSCALE = "-B${xanot}/${yanot}wSne";
       `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift -Y$yshift $BSCALE>>$PS_file`;
   }

   print "Now is $title[$ii] ...\n";
#   `psxy boundary.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W20/180/180/180 >> $PS_file`;

   `psxy $comp2[$ii] -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
   if($offset_S !=0.){
#        `psxy $comp3[$ii] -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Liner >> $PS_file`;
   }
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

#if($plot_raw == 1) {`rm @comp2 @comp3`;}

# ----------------- draw models -------------------------
if($draw_model == 1) {
   print "----- Plotting model -----\n";
   $XMIN = $x_o/$deg + 1.0;
   $XMAX = ($x_o + 1300)/$deg - 1;
   $YMIN = 0;
   $YMAX = 300;
   $velm = "vel.xy";
   $output_vd = 1; # if output of velocity model is of the difference between (c44-c55)/(c44+c55)/2

   &create_model("semmodel",2,"zone1.xy",$x_o,0);
   &model_layers("semmodel",$x_o);

   $yanot  = 50;
   $xanot  = 2;
 
   $XLEN   = 11.8;
   $YLEN   = -$XLEN/($XMAX-$XMIN)*($YMAX-$YMIN)/$deg*1.;
   $Line = "4/50/50/50";
 
   $PORT   = ""; 
   $SCALE  = "$XLEN/$YLEN";
   $RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
   $FSCALE = "-F255/255/255";
   $BSCALE = "-B${xanot}/${yanot}WSne";

   $B= $BSCALE;
   $Rslice="-R$XMIN/$XMAX/$YMIN/$YMAX";

   $CPT= "vel_gray.cpt";
   `cp ../$CPT $CPT`;
   $xshift = -5.8;
   $yshift = 18.;

    if($output_vd == 1){
      $cmax = 0.1;
      open(SVC,"|/home/zhl/work/tomo/software/x/svcpt13_table_cont_zhl");
      print SVC "-$cmax $cmax\n";
      close(SVC);
      `mv svel13.cpt $CPT`;
    }

   `surface $velm -Gvelm.grd -R$RSCALE -I0.01/0.1 -T0.25 -V`;
   `grdimage velm.grd -X$xshift -Y$yshift -JX$SCALE -R$RSCALE -C$CPT  -K -O >> $PS_file`;

#   `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  $BSCALE >>$PS_file`;
   `psxy layers.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W5/200/200/200  >> $PS_file`;

#   `psxy $int_model -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W5/200/0/0  >> $PS_file`;

   &PSTEXT($XMIN -($XMAX-$XMIN)/8 , $YMIN + ($YMAX - $YMIN) /2 ,11, 90, 1, 10, "Depth (km)",$PS_file);
   &PSTEXT($XMIN+($XMAX-$XMIN)/2 ,$YMAX + ($YMAX-$YMIN)/3,11, 0, 1, 10, "Distance (\260)",$PS_file);
   &PSTEXT($XMIN -1.3 , $YMIN - ($YMAX - $YMIN) /7 ,15, 0, 1, 10, "(a)",$PS_file);

   `psxy arrow.xy -JX$SCALE -R$RSCALE  -N -K -O $PORT -SV0.1c/0.2c/0.13c -G100/100/100 >> $PS_file`;

   open(RC,"receiver.xy");
   $nl=<RC>; chomp($nl);
   $longs=<RC>; chomp($longs);
   @segs=split(/ +/,$longs);
   open(SS, ">sta.xy");
   for($ii=0; $ii<$nl; $ii++){
       $x = ($segs[$ii] + $x_o)/$deg;
       print SS "$x -6\n";
   }
   close(SS);
   close(RC);

   `psxy sta.xy -JX$SCALE -R$RSCALE  -N -K -O $PORT -Si0.15 -G0/0/0 >> $PS_file`;
#   `psxy arrow.xy -JX$SCALE -R$RSCALE  -N -K -O $PORT -SV0.1c/0.2c/0.13c -G100/100/100 >> $PS_file`;
   `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  $BSCALE >>$PS_file`;
}

# -------------------------------------------------------
# draw split parameters
$plot_sws = 0;
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

sub create_model{
   my($fd_file,$layer,$an_file1,$x_o,$layn) = @_;
   #$layer: the structural layer number 
   $back = 1200;
   open(FF, "$fd_file");
   $l = <FF>; chomp($l);
   $l_num = $l;
   $lay = $layn + 1;
   $xstart = 0;
   $xstep = 10;


   open(AN,">$an_file1");
   print AN ">\n";

   for($ll=0;$ll<$l_num;$ll++){ 
      $l  = <FF>; chomp($l);
      $an = <FF>; chomp($an);
      $bn = <FF>; chomp($bn);
      
      if($ll == ($l_num -$lay-1)) {
         undef(@an_segs);
         undef(@bn_segs);
         @an_segs= split(/ +/,$an);
         @bn_segs= split(/ +/,$bn);
         $pp = @an_segs; $temp = $an_segs[$pp-2];
         print "line =$ll pp=$pp temp=$temp \n";         
         for($li=0;$li<@an_segs-1;$li++){
            if($an_segs[$li] <= $back){
                $lx = ($an_segs[$li]+$x_o)/111.2;
                $ly = $top - $bn_segs[$li];
                print AN "$lx $ly\n";
                if($li==0) {$lx0=$lx; $ly0=$ly;}
            }
         }

        $mmm = 0;
        $xxx = $xstart;
        while($xxx < $an_segs[$pp-2]){
            $tx[$mmm] = $xxx;
            for($li=0;$li<@an_segs-2;$li++){
               if($xxx>=$an_segs[$li] && $xxx<$an_segs[$li+1]){
                    $slope = ($bn_segs[$li+1]-$bn_segs[$li])/($an_segs[$li+1]-$an_segs[$li]);
                    $ty[$mmm] = $bn_segs[$li] + $slope*($xxx-$an_segs[$li]);
               }
            }
#            print "$mmm $xxx $ty[$mmm]\n";
            $mmm ++;
            $xxx = $xxx + $xstep;

        }
      }

      if($ll == ($l_num-$lay)) {
         print "last line =$lay\n";
         @stiffs = split(/ +/,$l);
         $rho = $stiffs[14];
         $vp = sqrt($stiffs[1]/$rho);
         $vs = sqrt($stiffs[10]/$rho);

         undef(@an_segs);
         undef(@bn_segs);

         @an_segs= split(/ +/,$an);
         @bn_segs= split(/ +/,$bn);
         
         for($li=0;$li<@an_segs;$li++){
            $point = @an_segs - 2 - $li;
            $lx = ($an_segs[$point]+$x_o)/$deg;
            $ly = $top - $bn_segs[$point];
            if($an_segs[$point] <= $back) { print AN "$lx $ly\n";}
         }

        for($kk=0;$kk<$mmm;$kk++){
            $xxx = $tx[$kk];
            for($li=0;$li<@an_segs-2;$li++){
               if($xxx>=$an_segs[$li] && $xxx<$an_segs[$li+1]){
                    $slope = ($bn_segs[$li+1]-$bn_segs[$li])/($an_segs[$li+1]-$an_segs[$li]);
                    $ty2[$kk] = $bn_segs[$li] + $slope*($xxx-$an_segs[$li]);
               }
            }
        }
         print AN "$lx0 $ly0\n";
         close(AN);
      }

   } 
   print AN "$lx0 $ly0\n";
   close(FF);
   close(AN);
   print "calculate offset vp=$vp vs=$vs\n";
   open(SLO,">offset");
   for($kk=0;$kk<$mmm;$kk++){
       $xxx = ($tx[$kk] + $x_o)/$deg;
       $time = ($ty2[$kk]-$ty[$kk])*(-1./$vp + 1./$vs);
       print SLO "$time $xxx \n";
#       print "$kk x=$tx[$kk] y=$ty[$kk] y1=$ty2[$kk] $xxx $time\n";
   }
   close(SLO);
}

sub model_layers{
   my($fd_file,$x_o) = @_;
   #$layer: the structural layer number 
   open(FF, "$fd_file");
   $l = <FF>; chomp($l);
   $l_num = $l;

   open(AN,">layers.xy");
   for($ll=0;$ll<$l_num;$ll++){
#   for($ll=$l_num-6;$ll<$l_num;$ll++){  
      $l  = <FF>; chomp($l);
      $an = <FF>; chomp($an);
      $bn = <FF>; chomp($bn);
      undef(@an_segs);
      undef(@bn_segs);

      @an_segs= split(/ +/,$an);
      @bn_segs= split(/ +/,$bn);
      print AN ">\n";
      for($npp=0;$npp<@an_segs;$npp++){
          $ax = ($an_segs[$npp]+$x_o)/$deg;
          $bx = $top - $bn_segs[$npp];
          print AN "$ax $bx\n";
      }  
   } 
   close(FF);
   close(AN);
}


