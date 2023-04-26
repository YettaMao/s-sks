# ---------------------------------------------------
#
# perl script for plotting U W components from GRT calculations
# 
# Written by Liang Zhao, IGG.CAS
# Jan 24, 2008
# ---------------------------------------------------
# see src/aserpsv_sem.c for details 

@ARGV == 0 || die "Usage: perl Plot_GRT_sh.pl ";
$input = "GRT_Plot_sh.xy";
$grtout_v = "grt_v.xy";
@comp  = ("V");
@title = ("v");
$one_degree = 111.2;
$enlarge = 1;

$PS_file = "grt_sh.ps";

$XMIN    = 0;
$XMAX    = 100;

## =====================================================
# read the output from asersh_sem.c -- GRT v
# & creat corrdinate file for GMT plot
$YMIN = 999;
$YMAX = -999;
$ampmax_v = -999;
open(GRTIN, $input);
$l = <GRTIN>; chomp($l);
@seg  = split(/ +/,$l);
$t_num = $seg[0]; $xtrace = $seg[1];
print "nt=$t_num xtrace= $xtrace \n";

for($ii=0;$ii<$xtrace;$ii++)
{
    $l = <GRTIN>; chomp($l);
    @seg  = split(/ +/,$l);
    # $is_v = 1: u; if =0: w
    $is_v = $seg[1];
#    print "$seg[0] $seg[1] $seg[2] is_v=$is_v\n";
    $dist = $seg[2]/$one_degree;
    if($dist < $YMIN) { $YMIN = $dist;}
    if($dist > $YMAX) { $YMAX = $dist;}
    for($jj=0; $jj<$t_num; $jj++)
    {
       $l = <GRTIN>; chomp($l);
       @seg = split(/ +/, $l);
       $time = $seg[0];
       if($ampmax_v < (abs($seg[1]))) {$ampmax_v = (abs($seg[1]));}
       if($ii==0 && $jj == 0) {$XMIN = $time;}
       if($ii==0 && $jj == $t_num-1) {$XMAX = $time;}
     }
}
close(GRTIN);

$amp_max = $ampmax_v;
print "maximum amplitude= $amp_max\n";

$Height = ($YMAX - $YMIN)/$xtrace;
$scale_u = $Height/$ampmax_v*$enlarge;
$scale   = $Height/$amp_max*$enlarge;
print "scale_v= $scale_u \n";
open(GRTIN, $input);
open(GRTU,"> $grtout_v");

$l = <GRTIN>; chomp($l);
for($ii=0;$ii<$xtrace;$ii++)
{
    $l = <GRTIN>; chomp($l);
    @seg  = split(/ +/,$l);
    $is_v = $seg[1];
    $dist = $seg[2]/$one_degree;
    if($is_v == 1) {
        print GRTU ">> \n";
    }
# filter here ------------------------------------

# ------------------------------------------------
    for($jj=0; $jj<$t_num; $jj++)
    {
       $l = <GRTIN>; chomp($l);
       @seg = split(/ +/, $l);
       $time = $seg[0];
       if($is_v == 1) {
           $uv = $seg[1] * $scale + $dist;
           print GRTU "$time $uv\n";
       }
    }
}
close(GRTIN);
close(GRTU);
close(GRTW);

## =======================================================

$YMIN -= 0.15; 
$YMAX += 0.15; 

print "xmin= $XMIN, xmax= $XMAX ; ymin= $YMIN, ymax= $YMAX \n";
$xshift = 2;
$yshift  = 5;
$yanot   = 1;
$xanot   = 50;

$YLEN       = "-20c";
$XLEN       = "7.5c";
$offset = "8c";

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-BWSne";
$BSCALE2 = "-BwSne";
$Line = "0/0/220";

$textsize = 15;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/20;
$y_offset = ($YMAX - $YMIN)/20;
print "plotting , please wait ...\n";
`gmt begin grt_sh ps`;
for($ii =0; $ii< @comp; $ii++)
{
   if($ii == 0) {`gmt basemap -JX$SCALE -R$RSCALE -X$xshift -Y$yshift $BSCALE  -Bxa${xanot} -Bya${yanot}`;}
   else { `gmt basemap -JX$SCALE -R$RSCALE -X$offset $BSCALE2 -Bxa${xanot} -Bya${yanot}`;}
 
   print "Now is $title[$ii] ...\n";

   if($ii == 0) {
      `gmt plot $grtout_v -JX$SCALE -R$RSCALE  -N -W1p,$Line`;
       $y = $YMIN - $y_offset;
       &PSTEXT($XMIN+($XMAX-$XMIN)/2,$y ,$textsize, $textangle, 1, 10, $title[$ii],$PS_file);
   }
}
`gmt end `;

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| gmt text -R$RSCALE -N -JX$SCALE -F+f+a+j");
        print GMT "$xx $yy $textsize,$textfont $textangle $just $text\n";
    close(GMT);
}
