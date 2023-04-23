# perl script for checking the validation of Lengenfre by Legendre() in Lagrange.h
# input: legendre.xz

$input = "/home/zhl/work/sem/applications/epi105/legendre.xz";
$PS_file = "legendre.ps";

$xshift = 4;
$yshift = 20;
$yanot  = 0.5;
$xanot  = 0.5;

$YLEN   = 6;
$XLEN   = 6;

$XMIN = -1;
$YMIN = -0.5;
$XMAX = 1;
$YMAX = 1.1;


$PORT   = "-P"; 
$SCALE  = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line   = "3/80/80/80";

`gmtset ANOT_FONT_SIZE 10`;
`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
`psxy $input -JX -R -M -N -K -O -P -W2  >> $PS_file`;

&PSTEXT($XMIN-0.5,($YMIN+$YMAX)/2 ,10,90, 1, 12,"value",$PS_file);
&PSTEXT(($XMIN+$XMAX)/2,$YMIN - 0.3 ,10, 0, 1, 12,"x ",$PS_file);

`psbasemap -JX -R $PORT -O -B >>$PS_file`;

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}
