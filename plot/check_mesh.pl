# perl script for checking the validation of meshing by discret() in sem2d.c
# input: grid.xy

$input = "/home/zhl/work/sem/applications/epi105/grid.xz";
$PS_file = "mesh.ps";

$xshift = 4;
$yshift = 20;
$yanot  = 100;
$xanot  = 200;

$YLEN   = -5;
$XLEN   = 12;

$XMIN = 0;
$YMIN = 0;
$XMAX = 1000;
$YMAX = 410;


$PORT   = "-P"; 
$SCALE  = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$Line   = "3/80/80/80";

`gmtset ANOT_FONT_SIZE 10`;
`surface $input -R$RSCALE -I0.5 -Gmesh.grd`;
`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;
`grdimage mesh.grd -R -JX -Cmesh.cpt -K -O >>$PS_file`;

`grdcontour  mesh.grd  -JX -Cmesh.cpt $PORT -A0.2+f10 -W8/255/255/255 -K -O >>$PS_file`;
&PSTEXT($XMAX+15,$YMAX-50 ,9,90, 1, 8,"Density (g/cm3)",$PS_file);

&PSTEXT($XMIN+330,$YMAX+100 ,12,0, 1, 8,"Distance (km)",$PS_file);
&PSTEXT($XMIN-160,$YMAX-20 ,12,90, 1, 8,"Depth (km)",$PS_file);

#`psscale -D12.6c/3c/4c/0.4c -L -N -Cmesh.cpt -I0.2 -B50 -O >>$PS_file`;

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}
