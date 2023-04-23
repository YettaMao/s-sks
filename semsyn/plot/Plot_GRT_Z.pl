# ---------------------------------------------------
#
# perl script for plotting U W components from GRT calculations
# 
# Written by Liang Zhao, IGG.CAS
# Jan 24, 2008
# ---------------------------------------------------
# see src/aserpsv_sem.c for details 

@ARGV == 0 || die "Usage: perl Plot_GRT.pl ";
$input = "GRT_Plot_z.xy";
$grtout_u = "grt_u_z.xy";
$grtout_w = "grt_w_z.xy";
@comp  = ("U","W");
@title = ("u","w");

$enlarge = 1;

$PS_file = "grtz.ps";

$XMIN    = 0;
$XMAX    = 100;

## =====================================================
# read the output from aserpsvfd_zl.c -- GRT u and v
# & creat corrdinate file for GMT plot
$YMIN = 999;
$YMAX = -999;
$ampmax_u = 0;
$ampmax_w = 0;
open(GRTIN, $input);
$l = <GRTIN>; chomp($l);
@seg  = split(/ +/,$l);
$t_num = $seg[0]; $xtrace = $seg[1];
print "nt=$t_num xtrace= $xtrace \n";

for($ii=0;$ii<$xtrace*2;$ii++)
{
    $l = <GRTIN>; chomp($l);
    @seg  = split(/ +/,$l);
    # $is_u = 1: u; if =0: w
    $is_u = $seg[1];
#    print "$seg[0] $seg[1] $seg[2]\n";
    $dist = $seg[2];
    if($dist < $YMIN) { $YMIN = $dist;}
    if($dist > $YMAX) { $YMAX = $dist;}
    for($jj=0; $jj<$t_num; $jj++)
    {
       $l = <GRTIN>; chomp($l);
       @seg = split(/ +/, $l);
       $time = $seg[0];
       if($is_u == 1 && $ampmax_u < (abs $seg[1])) {$ampmax_u = (abs $seg[1]);}
       if($is_u == 0 && $ampmax_w < (abs $seg[1])) {$ampmax_w = (abs $seg[1]);}
       if($ii==0 && $jj == 0) {$XMIN = $time;}
       if($ii==0 && $jj == $t_num-1) {$XMAX = $time;}
     }
}
close(GRTIN);

if($ampmax_u>$ampmax_w) {$amp_max=$ampmax_u;}
else {$amp_max=$ampmax_w;}

print "maximum amplitude= $amp_max\n";

$Height = ($YMAX - $YMIN)/$xtrace;
$scale_u = $Height/$ampmax_u*$enlarge;
$scale_w = $Height/$ampmax_w*$enlarge;
$scale   = $Height/$amp_max*$enlarge;
print "scale_u= $scale_u scale_w= $scale_w\n";
open(GRTIN, $input);
open(GRTU,"> $grtout_u");
open(GRTW,"> $grtout_w");
$l = <GRTIN>; chomp($l);
for($ii=0;$ii<$xtrace*2;$ii++)
{
    $l = <GRTIN>; chomp($l);
    @seg  = split(/ +/,$l);
    $is_u = $seg[1];
    $depth = $seg[2];
    if($is_u == 1) {
        print GRTU ">> \n";
    }
    else{
        print GRTW ">> \n";
    }
# filter here ------------------------------------

# ------------------------------------------------
    for($jj=0; $jj<$t_num; $jj++)
    {
       $l = <GRTIN>; chomp($l);
       @seg = split(/ +/, $l);
       $time = $seg[0];
       if($is_u == 1) {
           $uv = $seg[1] * $scale + $depth;
           print GRTU "$time $uv\n";
       }
       else{
           $uv = $seg[1] * $scale + $depth;
           print GRTW "$time $uv\n";
       }
    }
}
close(GRTIN);
close(GRTU);
close(GRTW);

## =======================================================

$YMIN -= 10; 
$YMAX += 10; 

print "xmin= $XMIN, xmax= $XMAX ; ymin= $YMIN, ymax= $YMAX \n";
$xshift = 2;
$yshift  = 5;
$yanot   = 100;
$xanot   = 50;

$YLEN       = -20;
$XLEN       = 7.5;
$offset = $XLEN + 0.5;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";
$BSCALE2 = "-B${xanot}/${yanot}wSne";
$Line = "3/0/0/220";

$textsize = 15;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/20;
$y_offset = ($YMAX - $YMIN)/20;
print "plotting , please wait ...\n";
for($ii =0; $ii< @comp; $ii++)
{
   if($ii == 0) {`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;}
   else { `psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$offset $BSCALE2>>$PS_file`;}
 
   print "Now is $title[$ii] ...\n";

   if($ii == 0) {
      `psxy $grtout_u -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
       $y = $YMIN - $y_offset;
       &PSTEXT($XMIN+($XMAX-$XMIN)/2,$y ,$textsize, $textangle, 1, 10, $title[$ii],$PS_file);
   }
   else {
      `psxy $grtout_w -JX$SCALE -M -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
       $y = $YMIN - $y_offset;
       &PSTEXT($XMIN+($XMAX-$XMIN)/2,$y ,$textsize, $textangle, 1, 10, $title[$ii],$PS_file);
   }
}

`psbasemap -JX -R -O -B>>$PS_file`;

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}
