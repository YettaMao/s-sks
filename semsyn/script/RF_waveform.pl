# perl script for plotting RF 
# plotting Z and R components
# Written by Liang Zhao

@ARGV == 2 || die "Usage: $0 RF PS_FILE\n";
my($rf, $PS_file) = @ARGV;

$Taup="/usr/local/TauP-2.1.1/bin/taup_time";
$phase= "P";
$tout="tout";
$pi = 3.1415926535897932;
$deg = 6371*$pi/180.;

$suffix="*.cz";

$XMIN    = -3;
$XMAX    = 17;
$X_O = -$XMIN;

`gmtset ANOT_FONT_SIZE 10`;
`gmtset ANOT_OFFSET 0.04`;
`gmtset FRAME_WIDTH 0.03`;

# change ZNE system to ZRT system
$here = `pwd`; chomp($here);

@Z_FILES = `ls $suffix`;
$SACFILE = $Z_FILES[0];
chomp($SACFILE);

print "the sacfile is $SACFILE \n";
$sdepth = `saclst evdp <$SACFILE`; chomp($sdepth);
if($sdepth>600) {$sdepth /=1000;}

$gcarc = `saclst gcarc < $SACFILE`; chomp($gcarc);
$baz  = `saclst baz < $SACFILE`; chomp($baz);

$min_g = 180;
$max_g = 0.;
for($i=0; $i<@Z_FILES; $i++)
{
   $SACFILE = $Z_FILES[$i];chomp($SACFILE);
   $gcarc[$i] = `saclst gcarc < $SACFILE`; chomp($gcarc[$i]);
   if($gcarc[$i] > $max_g) {$max_g = $gcarc[$i];}
   if($gcarc[$i] < $min_g) {$min_g = $gcarc[$i];}   
}

$g_span = $max_g -$min_g;

$YMIN    = $min_g - $g_span/30.;
$YMAX    = $max_g + $g_span/30.;
$YR= $YMAX-$YMIN;
if($YR < 0.01) {$YR=1;}

$xshift  = 1.1;
$yshift  = 2;
$yanot   = 1;
$xanot   = 5;

$YLEN    = -8.0;
$XLEN    = 2.6;
$unit    = $YR/2.56/$YLEN*3;

$PORT   = "-P"; 
$SCALE  = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE = "-F255/255/255";
$BSCALE = "-B${xanot}/${yanot}WSne";

`gmtset MEASURE_UNIT inch`;
# ---- plot radial components first -----
$textsize  = 13;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/12;
$y_offset = ($YMAX - $YMIN)/12;
print "\n\n plot radial component first, please wait ...\n";
`psbasemap -JX$SCALE -R$RSCALE -K $PORT  -X$xshift -Y$yshift $BSCALE>$PS_file`;

$textsize = 10;
$y = $YMAX + $y_offset/2.4;
&PSTEXT($XMAX+0.15*($XMAX+$XMIN)/2,$y ,12, 0, 1, 10, "time (sec)",$PS_file);
$y = ($YMAX + $YMIN)/2 + ($YMAX - $YMIN)/20;;
&PSTEXT($XMIN -($XMAX+$XMIN)/2.3,$y ,12, 90, 1, 10, "distance (degree)",$PS_file);
&PSTEXT($XMIN +($XMAX-$XMIN)/20,$YMIN+($YMAX-$YMIN)/50,9, 0, 3, 5, "radial",$PS_file);

# event information

$sac = $Z_FILES[0];chomp($sacz);
$evdp = `saclst evdp < $sac`; chomp($evdp);
if($evdp > 600) {$evdp/=1000;}
$evla = `saclst evla < $sac`; chomp($evla);
$evlo = `saclst evlo < $sac`; chomp($evlo);

$textsize = 8;
for($i=0; $i<@Z_FILES; $i++)
{
       $Line = "5/60/60/60";
       $sacz = $Z_FILES[$i];chomp($sacz);

       $stla = `saclst stla < $sacz`; chomp($stla);
       $stlo = `saclst stlo < $sacz`; chomp($stlo);
       `$Taup -ph $phase -h $evdp -evt $evla $evlo -sta $stla $stlo >$tout`;
       $ray_para= &read_TauP_time($tout)/$deg;

       if($i==0) {
           print "$Taup -ph $phase -h $evdp -evt $evla $evlo -sta $stla $stlo\n";      
           print "ray para= $ray_para\n";
       }

       @segs= split(/\./,$sacz);
       $sta = $segs[0];
       $ll=@segs;
       $segs[$ll-1]= "cn";
       $sacn= join(".",@segs);
       $segs[$ll-1]= "ce";
       $sace= join(".",@segs);
       $gcarc = `saclst gcarc < $sacz`; chomp($gcarc);


       open(SAC, "|sac");
       print SAC "r $sacz $sacn $sace\n";
       print SAC "ch a $amarker\n";
       print SAC "w over\n";
       print SAC "quit\n";
       close(SAC);

       $pp = $ray_para;
       print "$rf $sacz $sacn $sace $pp\n";
       `$rf $sacz $sacn $sace $pp`;
       $radial= $sta.".rad";
       $transverse= $sta.".tra";
       `mv radial.rf $radial`;
       `mv transverse.rf $transverse`;
 #      &PlotRF($radial,$gcarc,$unit);
 #     `psxy tmp -JX$SCALE -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
       &PlotRF_wiggle($radial,$gcarc,$unit);
       `pswiggle  tmp -JX -R -K -O $PORT -Z-0.01 -N -W$Line -G60 >> $PS_file`;
}

`rm tmp`;

# ---- plot transverse components next -----
print "\n plot transverse component , please wait ...\n";

$BSCALE = "-B${xanot}/${yanot}wSne";
$xshift = $XLEN*1.12;

`psbasemap -JX$SCALE -R$RSCALE -K -O $PORT  -X$xshift  $BSCALE>>$PS_file`;

&PSTEXT($XMIN +($XMAX-$XMIN)/20,$YMIN+($YMAX-$YMIN)/50,9, 0, 3, 5, "transverse",$PS_file);
$x = $ph2_time[$YMIN_int]; $y = $YMIN + $y_offset/8;
&PSTEXT($x,$y ,8, $textangle, 3, 5, $phase2,$PS_file);

for($i=0; $i<@Z_FILES; $i++)
{
       $sacz = $Z_FILES[$i];chomp($sacz);
       @segs= split(/\./,$sacz);
       $sta = $segs[0];

       $gcarc = `saclst gcarc < $sacz`; chomp($gcarc);
       $transverse= $sta.".tra";
       print "plotting $transverse\n";
 #      &PlotRF($transverse,$gcarc,$unit);
 #     `psxy tmp -JX$SCALE -R$RSCALE  -N -K -O $PORT -W$Line >> $PS_file`;
       &PlotRF_wiggle($transverse,$gcarc,$unit);
      `pswiggle  tmp -JX -R -K -O $PORT -Z-0.01 -N -W$Line -G60 >> $PS_file`;
       &PSTEXT($XMAX + $x_offset, $gcarc,8, 0, 1, 10, $sta,$PS_file);

}
`rm tmp`;
`psbasemap -JX$SCALE -R$RSCALE -O $PORT  $BSCALE>>$PS_file`;

##-----------------------------------------------------------
## subroutines following
sub PlotRF{
    my($file,$gcarc,$unit)=@_;
    open(RR, $file);
    $line=<RR>;
    $line=<RR>;
    $line=<RR>;
    open(TMP,">tmp");
    $num = 0;
    while($line=<RR>){
        chomp($line);
        @segs= split(/ +/,$line);
        $xx= $segs[0];
        $yy= $gcarc + $segs[1]*$unit;
        print TMP "$xx $yy\n";
#        if($num==0) {print "$xx $yy\n";}
        $num ++;
    }
    close(RR);
    close(TMP);
}

sub PlotRF_wiggle{
    my($file,$gcarc,$unit)=@_;
    open(RR, $file);
    $line=<RR>;
    $line=<RR>;
    $line=<RR>;
    open(TMP,">tmp");
    $num = 0;
    while($line=<RR>){
        chomp($line);
        @segs= split(/ +/,$line);
        $xx= $segs[0];
        $yy=  - $segs[1]/40;
        print TMP "$xx $gcarc $yy\n";
#        if($num==0) {print "$xx $yy\n";}
        $num ++;
    }
    close(RR);
    close(TMP);
}

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| pstext -R$RSCALE -K -O -N -JX$SCALE >> $ps");
        print GMT "$xx $yy $textsize $textangle $textfont $just $text\n";
    close(GMT);
}

sub read_TauP_time{
   my($in)=@_;
   open(IN,$in);
   $line= <IN>; chomp($line);
   $line= <IN>; chomp($line);
   $line= <IN>; chomp($line);
   $line= <IN>; chomp($line);
   $line= <IN>; chomp($line);
   $line= <IN>; chomp($line);
   @segs= split(/ +/,$line);
   $amarker = $segs[4];
   $ray_para = $segs[5];
#   print "3:$segs[3] 4=$segs[4] 5=$segs[5]\n";
   close(IN);
   return $ray_para;
}
# 
#Model: iasp91
#Distance   Depth   Phase   Travel    Ray Param  Takeoff  Incident  Purist    Purist
#  (deg)     (km)   Name    Time (s)  p (s/deg)   (deg)    (deg)   Distance   Name 
#-----------------------------------------------------------------------------------
#   70.32    10.0   P        673.72     6.123     18.65    18.62    70.32   = P    


