# ---------------------------------------------------
#
# perl script for plotting U V W components 
# 
# Written by Liang Zhao, IGG.CAS
# 20141115
# ---------------------------------------------------

$plot_w = 1;
$offset = 1;
$draw_model = 1;
$x_o = 9000;
$pi = 3.14159265;
$deg = 6371*$pi/180.;

$sks_t_offset = 0;

$column_num = 2; 
if($plot_w == 1) {
   @comp  = ("u.xy", "v.xy", "w.xy");
   @title = ("Radial","Transverse","Vertical");
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
#$YMIN = $seg[0]-0.2; 
#$YMAX = $seg[1]+0.2;
#$XMIN += 10;
#$XMAX -= 20;
$YMIN = $seg[0]-1; #20211103 maoyt
$YMAX = $seg[1]+0.5;
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
$BSCALE = "-BWSne";
$Line = "3/80/80/250";
$Line = "80/80/250";
$Time_line  = "4/0/0/0t20_20:15";
$STime_line = "4/0/0/255t20_20:15";
$Time_line  = "0/0/0";
$STime_line = "0/0/255";

$textsize = 15;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/20;
$y_offset = ($YMAX - $YMIN)/40;
print "plotting , please wait ...\n";
`gmt begin Fig1 ps`;
for($ii =0; $ii< @comp; $ii++)
{
   if($ii==0) {
         `gmt basemap -JX$SCALE -R$RSCALE -X$xshift -Y$yshift $BSCALE -Bxa$xanot -Bya$yanot`;
          &PSTEXT($XMAX + 40, $YMAX + $y_offset*2.5 ,10, 0, 1, 10, "Time (s)",$PS_file);
          &PSTEXT($XMIN - 6.5* $x_offset, $YMIN +($YMAX - $YMIN)/2,10, 90, 1, 10, "Distance (\260)",$PS_file);
       #   &PSTEXT($XMIN - 3.5* $x_offset, $YMIN -($YMAX - $YMIN)/20,17, 0, 1, 10, "(b)",$PS_file);#maoyt comment
   }
   else 
   {  
       $xshift = "7.3c"; $yshift = 0.;
       if($plot_w == 1) {
            $xshift = "5.3c"
       }

       $BSCALE = "-BwSne";
       `gmt basemap -JX$SCALE -R$RSCALE -X$xshift -Y$yshift $BSCALE -Bxa$xanot -Bya$yanot`;
   }
   print "Now is $title[$ii] ...\n";

   `gmt plot $comp[$ii] -W$Line`;
   if($offset <=0) {
       `gmt plot time.of -N -W$Time_line `;
       `gmt stime.of -N  -W$STime_line`;
   }

   $y = $YMIN + $y_offset*0.6;
   &PSTEXT($XMIN+($XMAX-$XMIN)/2,$y ,9, $textangle, 1, 10, $title[$ii],$PS_file);

   `gmt basemap $BSCALE`;
}


# ----------------- draw models -------------------------
if($draw_model == 1) {
   print "----- Plotting model -----\n";
   $XMIN = $x_o/$deg + 1.0;
   $XMAX = ($x_o + 1300)/$deg - 1;
   $YMIN = 0;
   $YMAX = 300;
   $velm = "vel.xy";
   $output_vd = 0; # if output of velocity model is of the difference between (c44-c55)/(c44+c55)/2

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
   $BSCALE = "-BWSne";

   $B= $BSCALE;
   $Rslice="-R$XMIN/$XMAX/$YMIN/$YMAX";

   $CPT= "vel_gray.cpt";
   `cp ../$CPT $CPT`;
   $xshift = -8.8;
   $yshift = 18.;

    if($output_vd == 1){
      $cmax = 0.1;
      open(SVC,"|/home/maoyt/work/tomo/software/x/svcpt13_table_cont_zhl");
      print SVC "-$cmax $cmax\n";
      close(SVC);
      `mv svel13.cpt $CPT`;
    }
   `gmt surface $velm -Gvelm.grd -R$RSCALE -I0.01/0.1 -T0.25 `;
   `gmt grdimage velm.grd -X$xshift -Y$yshift -JX$SCALE -R$RSCALE -C$CPT `;

   `gmt basemap -JX$SCALE -R$RSCALE $BSCALE -Bxa$xanot -Bya$yanot`;
   `gmt plot layers.xy -JX$SCALE -R$RSCALE -W200/200/200 `;

#   `gmt plot $int_model -JX$SCALE -R$RSCALE -W200/0/0 `;

   &PSTEXT($XMIN -($XMAX-$XMIN)/8 , $YMIN + ($YMAX - $YMIN) /2 ,11, 90, 1, 10, "Depth (km)",$PS_file);
   &PSTEXT($XMIN+($XMAX-$XMIN)/2 ,$YMAX + ($YMAX-$YMIN)/3,11, 0, 1, 10, "Distance (\260)",$PS_file);
   #&PSTEXT($XMIN -1.3 , $YMIN - ($YMAX - $YMIN) /7 ,15, 0, 1, 10, "(a)",$PS_file);

   #  `gmt plot arrow.xy -SV0.1c/0.2c/0.13c -G100/100/100`;

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
   `gmt plot sta.xy -N -Si0.15c -Gblack`;

#   `gmt plot arrow.xy -SV0.1c/0.2c/0.13c -G100/100/100 `;
}
`gmt end`;
#`gs $PS_file`;


sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| gmt text -R$RSCALE -N -F+f+a+j ");
        print GMT "$xx $yy $textsize,$textfont $textangle $just $text\n";
    close(GMT);
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

sub S_tab{
   my($gcarc) = @_;
   open(TT,"taup_curve_S.gmt");
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


