# perl script for checking the validation of GRT-SEM interface
# input: inter.dat

$input = "/home/zhl/work/sem/applications/epi90/inter.dat";
$input2 = "/home/zhl/work/sem/applications/epi90/inter2.dat";
$input3 = "/home/zhl/work/sem/applications/epi90/inter3.dat";
$PS_file = "interface.ps";

$ymax= 0;
$num = 0;
open(RR,"$input");
while($line= <RR>){
   chomp($line);
   @segs= split(/ +/, $line);
   $vals = $segs[1];
   if(fabs($vals) > $ymax) {$ymax = fabs($vals);}
   if($num == 0) {print "val 0= $vals \n";}
   $num ++ ;
}
close(RR);

print "ymax= $ymax\n";

$xshift = 3;
$yshift = 12;
$yanot  = $ymax / 2.;
$xanot  = 50;

$YLEN   = 6;
$XLEN   = 12;

$XMIN = 0;
$YMIN = -$ymax* 1.1;
$XMAX = $XMIN + 80;
$YMAX = $ymax * 1.1;


$PORT   = "-P"; 
$SCALE  = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line   = "5/80/80/80";

`gmtset ANOT_FONT_SIZE 10`;
`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
`psxy $input -JX -R -N -K -O -P -W2  >> $PS_file`;

`psbasemap -JX -R -K $PORT -O -X -Y8 $BSCALE>>$PS_file`;
`psxy $input2 -JX -R -N -K -O -P -W2/0/0/255  >> $PS_file`;


$YMIN = -$ymax* 0.1;
$YMAX = $ymax * 0.1;
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
`psbasemap -JX -R$RSCALE -K $PORT -O -X -Y-16 $BSCALE>>$PS_file`;
`psxy $input3 -JX -R -N -K -O -P -W2/0/255/255  >> $PS_file`;

&PSTEXT($XMIN-0.5,($YMIN+$YMAX)/2 ,10,90, 1, 12,"value",$PS_file);
&PSTEXT(($XMIN+$XMAX)/2,$YMIN - 0.3 ,10, 0, 1, 12,"x ",$PS_file);

`psbasemap -JX -R $PORT -O -B >>$PS_file`;

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}

sub fabs{
    my($value) = @_;
    $val = $value;
    if($val < 0) {$val = $val* -1;}
    return $val;
}
