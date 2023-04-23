# perl script for plot the output from aserpsv_fem

@comp  = ("grt.u", "grt.w");
$PS_file = "grt_b.ps";

$XMIN=0; $XMAX=150;
$YMIN=100; $YMAX=103.;

$xshift = 2.5;
$yshift = 15;
$yanot  = 1;
$xanot  = 30;

$YLEN   = -12;
$XLEN   = 6;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WS";
$Line = "1/30/30/30";
$Time_line = "4/0/0/0t20_20:15";
$mline = "4/100/180/180";

print "plotting , please wait ...\n";
for($ii =0; $ii< @comp; $ii++)
{
   if($ii==0) {
       `psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
        &PSTEXT($XMIN+($XMAX-$XMIN)/15, $YMIN + 0.2 ,10, 0, 1, 10, "x1",$PS_file);
   }
   else 
   {  
       $xshift2 = $XLEN + 0.8;
       $BSCALE = "-B${xanot}/${yanot}S";
       `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift2 -Y0 $BSCALE>>$PS_file`;
        &PSTEXT($XMIN+($XMAX-$XMIN)/15, $YMIN + 0.2 ,10, 0, 1, 10, "x2",$PS_file);
   }

   `psxy $comp[$ii] -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$mline >> $PS_file`;
}

`psbasemap -JX -R -O -B>>$PS_file`;

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}
