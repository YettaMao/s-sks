# perl scripts for drawing snapshot of fd results
# Written by Liang Zhao, IGG, CAS
# Aug 17, 2006

@ARGV == 4 || die "Usage: $0 filename Figer_index snap_time spatial_step\n";
my($filename,$figer_index,$snap_time,$h) = @ARGV;

$TITLE = "(".$figer_index.")";
$m_time=" t= ".$snap_time." s";

$maxv = -999;
$maxv2 = -999;
@comp  = ("w.xy", "u.xy");
$PS_file = "snapshot.ps";
$modelfile = "semrf";
$m_xyfile = "model.xy";
$enlarge = 30.0;
$skip = 3;

`gmtset MEASURE_UNIT cm`;

open(INPUT,"$filename");
$l=<INPUT>; chomp($l);
@seg = split(" ",$l);
$x[0] = $seg[0];
$y_min = $seg[1];
   $y[0] = $seg[1];
   $u[0] = $seg[2];
   $v[0] = $seg[3];
$new_x = 0;
$raw_num = 0;
$col_num = 0;
$num = 0;
while($l=<INPUT>){
   $y_old = $seg[1];
   $num ++;
   chomp($l);
   @seg = split(" ",$l);
   $y[$num] = $seg[1];
   $u[$num] = $seg[2];
   $v[$num] = $seg[3];
   if(&fabs($seg[2])<1e20){
      if($maxv < &fabs($seg[2])) {$maxv = &fabs($seg[2]);}
      if($maxv < &fabs($seg[3])) {$maxv = &fabs($seg[3]);}
#      if($maxv2 < &fabs($seg[3])) {$maxv2 = &fabs($seg[3]);}
      if($seg[0] > $x[$col_num]) {
         if($new_x == 0) {
             $raw_num = $num ;
             $y_max = $y_old;
         }
         $new_x = 1;
         $col_num ++;
         $x[$col_num] = $seg[0];
      }
   }
}
close(INPUT);
$col_num ++;

if($maxv ==0 ) {$maxv = 0.00001;}
if($maxv2 ==0 ) {$maxv2 = 0.00001;}

print STDOUT"raw=$raw_num col=$col_num\n";
print "x[0]= $x[0] x_max= $x[$col_num-1] ymin= $y_min ymax= $y_max\n";
print "maxv= $maxv maxv2=$maxv2\n";
$XMIN = $x[0]*$h; $XMAX = $x[$col_num-1]*$h;
$YMIN = $y_min*$h;$YMAX = $y_max*$h;   

$YMIN = -1;

$unit = ($XMAX-$XMIN)/$col_num;

# -----------------------------------------------------------------------
$xshift = 2.5;
$yshift = 20;
$yanot  = 1000;
$xanot  = 3000;

$YLEN   = 7/($XMAX-$XMIN)*($YMAX-$YMIN);
$XLEN   = 7;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}s";
$Line = "1/30/30/30";
$Time_line = "4/0/0/0t20_20:15";
$mline = "8/180/180/180";

print "plotting , please wait ...\n";
print "creating model\n";
&fdmodel($modelfile, $m_xyfile);
print "plotting snapshot\n";
for($ii =0; $ii< @comp; $ii++)
{
   if($ii==0) {
       `psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
        &PSTEXT($XMIN+($XMAX-$XMIN)/15, $YMIN + 100 ,10, 0, 1, 10, "x1",$PS_file);
        &PSTEXT($XMIN + 50, $YMAX + 190 ,16, 0, 1, 10, $TITLE,$PS_file);
        &PSTEXT($XMIN + 450, $YMAX + 140 ,10, 0, 1, 10, $m_time,$PS_file);
   }
   else 
   {  
       $xshift2 = $XLEN + 0.56;
       $BSCALE = "-B${xanot}/${yanot}s";
       `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift2 -Y0 $BSCALE>>$PS_file`;
        &PSTEXT($XMIN+($XMAX-$XMIN)/15, $YMIN + 100 ,10, 0, 1, 10, "x2",$PS_file);
   }

   open(TT, ">temp.xy");
   for($cl=0;$cl<$col_num;$cl++){
      $is_skip = &mod($x[$cl], $skip);
      if(! $is_skip) {
 #         print "column No. $cl $is_skip\n";
          print TT ">\n";
          for($rr=0;$rr<$raw_num;$rr++){
               $offset = $cl*$raw_num+$rr;
               $y_xy = $y[$offset]* $h;
               if($ii==0) {$amplitude = $x[$cl]*$h + $u[$offset]/$maxv*$unit*$enlarge;}
               if($ii==1) {$amplitude = $x[$cl]*$h + $v[$offset]/$maxv*$unit*$enlarge;}
               print TT "$amplitude $y_xy\n";
          }
      }
   }
   close(TT);
   `psxy $m_xyfile -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$mline >> $PS_file`;
   `psxy temp.xy -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
}

`psbasemap -JX -R -O -B>>$PS_file`;

`rm temp.xy`;
`rm $m_xyfile`;
#-------------------------------------------
sub fdmodel{
   my($modelfile, $m_xyfile) = @_;
   print "$modelfile, $m_xyfile\n";
   open(MM, "$modelfile");
   open(XY, ">$m_xyfile");
   $line = <MM>; chomp($line);
   @seg = split(" ", $line);
   $layer = $seg[0];
   print "nl= $layer\n";
   $iso = 1;
   for($ml=0; $ml< $layer; $ml++){
#      print "ml = $ml";
      $line = <MM>; chomp($line);
      @seg = split(" ", $line);
      $np= $seg[0]; $aniso = $seg[3];
      if($aniso != 0 && $iso == 1) {
          print XY ">\n";
          $iso = 0;
          for($nn=0; $nn< @an; $nn++){
             print XY "$an[$nn] $bn[$nn]\n";
          }
      }
      $line = <MM>; chomp($line);
      @an = split(" ", $line);
      $line = <MM>; chomp($line);
      @bn = split(" ", $line);
      if($aniso != 0) {
          print XY ">\n";
          for($nn=0; $nn< $np; $nn++){
             print XY "$an[$nn] $bn[$nn]\n";
          }
      }
   }
   close(MM);
   close(XY);
}

sub fabs{
   my($value) = @_;
   if($value > 0) {$abs = $value; }
   else {$abs = - $value;}
   return $abs;
}

sub mod{
   my($value, $div) = @_;
   $tempv = $value;
   while($tempv >= $div) {
      $tempv -= $div;
   }
   return $tempv;
}

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}
