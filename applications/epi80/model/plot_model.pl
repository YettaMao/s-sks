# perl scripts for plotting models and splitting results

@ARGV == 3 || die "Usage: $0 model_file width depth\n";
my($model_file,$width,$depth) = @ARGV;

$PS_file = "model.ps";
$temp = "temp_model";
$station = "../receiver.xy";
$phi_0 = 250;
$phi_180 = 100;
$BNS = 104.8; #degree

$xo = 11076;

$xshift = 2.5;
$yshift = 20;
$yanot  = 100;
$xanot  = 1;

$xbeg = 350;
$XMIN = ($xbeg+$xo)/111.2; $XMAX = ($xo+$width)/111.2;
$YMIN = 0; $YMAX = 220;
$YLEN   = -12*($YMAX-$YMIN)/($width-$xbeg);
$XLEN   = 12;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line = "3/30/30/30";
$Lines[0] = "5/30/30/30";
$Lines[1] = "5/50/50/50";
$Lines[2] = "5/70/70/70";
$Lines[3] = "5/90/90/90";
$Lines[4] = "5/110/110/110";
$Lines[5] = "5/120/120/120";

print "plotting , please wait ...\n";
print "creating model\n";
print "$XMIN $XMAX $YMIN $YMAX\n";

`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
&PSTEXT($XMIN+($XMAX-$XMIN)/2,$YMAX + ($YMAX - $YMIN)/7.8 ,10, 0, 1, 10, "Distance (\260)",$PS_file);
&PSTEXT($XMIN-($XMAX-$XMIN)/8,$YMIN +($YMAX - $YMIN)/2 ,10, 90, 1, 10, "Depth (km)",$PS_file);

&PSTEXT($BNS,$YMIN-15,8, 90, 1, 10, "BNS",$PS_file);

# Plot station 
&PlotSta($station,-3);
`psxy temp.xy  -JX$SCALE  -R$RSCALE  -N -K -O $PORT -W$Line -Si0.1  >> $PS_file`;

# Plot BNS
$BNS_line = "10/200/200/200";
&PlotBNS($BNS);
`psxy temp.xy  -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$BNS_line  >> $PS_file`;

open(MM, "$model_file");
$ll = <MM>; chomp($ll);
$lay_k =0;
for($ii=0;$ii<$ll;$ii++){
    $mark= <MM>; chomp($mark);
    $lan= <MM>; chomp($lan);
    $lbn= <MM>; chomp($lbn);
    @an = split(/ +/,$lan);
    @bn = split(/ +/,$lbn);
    if($mark=~/phi/){
#       print "$mark\n";
       $lay_k ++;
       @segs = split(/ +/,$mark);
       $phi[$ii] = $segs[2];
       $ggg = $phi_0 - $phi[$ii]*($phi_0-$phi_180)/180;
       @gsegs = split(/\./,$ggg);
       $grey = $gsegs[0];
       print "phi = $phi[$ii] grey=$grey\n";
       @last_an = split(/ +/,$ll_an);
       @last_bn = split(/ +/,$ll_bn);
       $num = 0;
       for($x=0;$x<$width;$x++){
          for($jj=1;$jj<@an;$jj++){
             if($x > $an[$jj-1] && $x < $an[$jj]){
                  $slope = ($bn[$jj]-$bn[$jj-1])/($an[$jj] - $an[$jj-1]);
                  $y1 = $bn[$jj-1] + $slope*($x-$an[$jj-1]);
             }
          }
          for($jj=1;$jj<@last_an;$jj++){
             if($x > $last_an[$jj-1] && $x < $last_an[$jj]){
                  $slope = ($last_bn[$jj]-$last_bn[$jj-1])/($last_an[$jj] - $last_an[$jj-1]);
                  $y2 = $last_bn[$jj-1] + $slope*($x-$last_an[$jj-1]);
             }
          }
          if($y1 > $y2){
             if($num == 1) {print "$x $y1 $y2\n";}
             $y_1[$num] = $y1;
             $y_2[$num] = $y2;
             $x_v[$num] = $x;
             $num ++;
          }
       }
       open(TT, ">$temp");
       print TT ">\n";
if($lay_k >1){
       for($jj=0;$jj<$num;$jj++){
           $xp = ($x_v[$jj]+$xo)/111.2;
           $yp = $depth - $y_2[$jj];
           print TT "$xp $yp\n";
       }
}
       for($jj=0;$jj<$num;$jj++){
           $kk = $num - $jj -1;
           $xp = ($x_v[$kk]+$xo)/111.2;
           $yp = $depth - $y_1[$kk];
           print TT "$xp $yp\n";
       }
       $xp = ($x_v[0]+$xo)/111.2; 
       $yp = $depth - $y_2[0];      
       print TT "$xp $yp\n";
       close(TT);
    }
    else {$phi[$ii] = 0.;}
    
    $fill = "$grey/$grey/$grey";
    `psxy $temp -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Lines[$ll]  >> $PS_file`;
    $ll_an = $lan;
    $ll_bn = $lbn
}
close(MM);

sub PlotSta{
    my($sta_file,$yy) =@_;
    open(ST,"$sta_file");
    open(REC,">temp.xy");
    $line = <ST>;
    $line = <ST>; chomp($line);
    @segs = split(/ +/,$line);
    for($rr=0; $rr<@segs;$rr++){
        print REC "$segs[$rr] $yy\n";
#        print "$segs[$rr] $yy\n";
    }
    close(REC);
    close(ST);
}

sub PlotBNS{
   my($xp) = @_;
   open(BB, ">temp.xy");
   print BB ">\n";
   print BB "$xp $YMAX\n";
   print BB "$xp $YMIN\n";
   close(BB);
}

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}
