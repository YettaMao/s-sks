# perl script for testing waveform.c using a sac file
# plotting R and T components

@ARGV == 3 || die "Usage: $0 DIR WF_C PS_FILE\n";
my($DIR, $WF, $PS_file) = @ARGV;
$ne2rt = "/home/maoyt/work/aniso/sws/local/PlotPl/NE2RT.pl";
$WFT   = "/home/maoyt/work/aniso/sws/local/bin/waveform_scale";

$temp = "temp.txt";

$caltime = "yes";
$phase = "S";
$align_phase = "SKSac";
$phase2 = "S"; # the 2nd phase for consideration

$XMIN    = -15;
$XMAX    = 45;
$X_O = -$XMIN;

$Old_dir = `pwd`; chomp($Old_dir);

`gmtset ANOT_FONT_SIZE 10`;
`gmtset ANOT_OFFSET 0.04`;
`gmtset FRAME_WIDTH 0.03`;

# change ZNE system to ZRT system
chdir $DIR;
$here = `pwd`; chomp($here);
print "Jump from $Old_dir -> $here\n";
`perl $ne2rt`;
chdir $Old_dir;

@R_FILES = `ls $DIR/*.BHR`;
@T_FILES = `ls $DIR/*.BHT`;
$SACFILE = $R_FILES[0];
chomp($SACFILE);

print "the sacfile is $SACFILE \n";
my($SACFILE,$sdepth,$nzyear,$nzjday,$nzhour,$nzmin,$nzsec) = split/\s+/,`saclst evdp nzyear nzjday nzhour nzmin nzsec f $SACFILE`;
if($sdepth>600) {$sdepth /=1000;}

my($SACFILE,$gcarc,$baz) = split/\s+/,`saclst gcarc baz f $SACFILE`; 
print "$nzyear $nzjday $nzhour  \n";

$min_g = 180;
$max_g = 0.;
for($i=0; $i<@R_FILES; $i++)
{
   $SACFILE = $R_FILES[$i];chomp($SACFILE);
   my($SACFILE,$gca) = split/\s+/,`saclst gcarc f $SACFILE`; 
   $gcarc[$i]=$gca;
   if($gcarc[$i] > $max_g) {$max_g = $gcarc[$i];}
   if($gcarc[$i] < $min_g) {$min_g = $gcarc[$i];}   
}
for($i=0; $i<@T_FILES; $i++)
{
   $SACFILE = $T_FILES[$i];chomp($SACFILE);
   my($SACFILE,$gca_n) = split/\s+/,`saclst gcarc f $SACFILE`; 
   $gcarc_n[$i]=$gca_n;
   if($gcarc_n[$i] > $max_g) {$max_g = $gcarc_n[$i];}
   if($gcarc_n[$i] < $min_g) {$min_g = $gcarc_n[$i];}   
}
$g_span = $max_g -$min_g;

$YMIN    = $min_g - $g_span/20.;
$YMAX    = $max_g + $g_span/20.;

for($i=0;$i<180;$i++){
    if($i<$YMIN && $i + 1 > $YMIN) {$YMIN_int = $i;}
}
open(TL,">time.of");
print TL "0. $YMIN\n";
print TL "0. $YMAX\n";
close(TL);

$xshift  = 3.1;
$yshift  = 2;
$yanot   = 1;
$xanot   = 15;

$YLEN       = -15.0;
$XLEN       = 6;

$PORT       = "-P"; 
$SCALE      = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE     = "-F255/255/255";
$BSCALE = "-BWSne";

# calculate time curves
if($caltime =~ /yes/)
{
    &caltimetable($phase,$sdepth);
    print INFO "cal time finished \n";
    # use hash %tcurves
    %tcurves = ();
    &timecurves(0,0);
    $dist = $tcurves{"dist"};
    while(($ph,$ts) = each %tcurves)
    {
	if($ph =~ /^$align_phase$/)
	{
	    $time = $ts;
        }
    }    
    @dis   = split(",", $dist);
    @times = split(",",$time);
    # create arrival time for 2nd phases related to align-phase
    open(TTT,">$phase2.of");
    while(($ph,$ts) = each %tcurves) {
	if($ph =~ /^$phase2$/){
	    $time_ph2 = $ts;
        }
    }    
    @times_ph2 = split(",",$time_ph2);
    for($j=0;$j<@dis;$j++){
         if($times_ph2[$j] > 0) {
            $ph2_time[$j] = $times_ph2[$j] - $times[$j];
            printf TTT "$ph2_time[$j] $j\n";
         }
         else {printf TTT ">>\n";}
    }
    close(TTT);
    # for R component
    for($i=0; $i<@gcarc; $i++)
    {
        $find = 0;
        for($j=2;$j<@dis; $j++)
        {
           if($gcarc[$i]<$dis[$j] && $gcarc[$i]> $dis[$j -1] && $find ==0)
           {
               $find = 1;
               $slope = ($times[$j]-$times[$j-1])/($dis[$j] - $dis[$j-1]);
               $a_times[$i] = $times[$j-1] + ($gcarc[$i]-$dis[$j -1])*$slope;
           }
        }
    }
    # for T component
    for($i=0; $i<@gcarc_n; $i++)
    {
        $find = 0;
        for($j=2;$j<@dis; $j++)
        {
           if($gcarc_n[$i]<$dis[$j] && $gcarc_n[$i]> $dis[$j -1] && $find ==0)
           {
               $find = 1;
               $slope = ($times[$j]-$times[$j-1])/($dis[$j] - $dis[$j-1]);
               $a_times_n[$i] = $times[$j-1] + ($gcarc_n[$i]-$dis[$j -1])*$slope;
           }
        }
    }    
}

`gmtset MEASURE_UNIT inch`;
#`gmt begin $PS_file ps`;
`gmt begin sks_test ps`;
# ---- plot R components first -----
$textsize  = 13;
$textangle = 0;
$x_offset = ($XMAX - $XMIN)/12;
$y_offset = ($YMAX - $YMIN)/12;
print "plot R component first, please wait ...\n";
`gmt basemap -JX$SCALE -R$RSCALE -X$xshift -Y$yshift $BSCALE -Bxa$xanot -Bya$yanot`;
`gmt plot time.of -JX$SCALE -R$RSCALE  -N  -W90/90/90 `;
`gmt plot $phase2.of -JX$SCALE -R$RSCALE  -N -W100/100/100 `;

$y = $YMIN - $y_offset;
if($nzjday<10) {$nzjday = "00".$nzjday;}
else {
   if($nzjday<100) {$nzjday = "0".$nzjday;}
}
if($nzhour<10) {$nzhour = "0".$nzhour;}
if($nzmin<10)  {$nzmin  = "0".$nzmin;}
if($nzsec<10)  {$nzsec  = "0".$nzsec;}
$title = "Event  ".$nzyear.":".$nzjday.":".$nzhour.":".$nzmin.":".$nzsec;
$info = "Event depth = ".$sdepth." km,   "."baz ~ ".$baz." degree";
&PSTEXT($XMAX+0.22*($XMAX+$XMIN)/2,$y ,$textsize, $textangle, 1, 10, $title,$PS_file);
$y = $YMIN - $y_offset/2.5;
$textsize = 10;
&PSTEXT($XMAX+0.16*($XMAX+$XMIN)/2,$y ,$textsize, $textangle, 1, 10, $info,$PS_file);
$y = $YMAX + $y_offset/2.4;
&PSTEXT($XMAX+0.15*($XMAX+$XMIN)/2,$y ,12, 0, 1, 10, "time (sec)",$PS_file);
$y = ($YMAX + $YMIN)/2 + ($YMAX - $YMIN)/20;;
&PSTEXT($XMIN -($XMAX+$XMIN)/2.3,$y ,12, 90, 1, 10, "distance (degree)",$PS_file);
&PSTEXT($XMIN +($XMAX-$XMIN)/20,$YMIN+($YMAX-$YMIN)/50,9, 0, 3, 5, "Radial",$PS_file);
$x = $ph2_time[$YMIN_int]; $y = $YMIN + $y_offset/8;
&PSTEXT($x,$y ,8, $textangle, 3, 5, $phase2,$PS_file);

$textsize = 8;
for($i=0; $i<@R_FILES; $i++)
{
    $Line = "0/0/0";
    $SACFILE = $R_FILES[$i];chomp($SACFILE);
    @tsegs = split(/_/,$SACFILE); 
    $sta_name = substr($tsegs[1],0,3);
    $B_time = $a_times[$i] - $X_O;
    $E_time = $B_time + ($XMAX -$XMIN);
    @name = split(/\./,$SACFILE);
    my($SACFILE,$long) = split/\s+/,`saclst kstnm f $SACFILE`; 
    @sta_name = split(/\./,$long);
    system("$WF $SACFILE $B_time $E_time $XMIN $YMAX $XMAX $YMIN 1 $temp ");
    open(TT,"scale.tmp");
    $unit_r[$i] = <TT>; chomp($unit_r[$i]);
#    print "No.$i scale = $unit_r[$i]\n";
    close(TT);
    `gmt plot $temp -JX$SCALE -R$RSCALE  -N -W$Line `;
#    &PSTEXT($XMAX + $x_offset, $gcarc[$i],$textsize, $textangle, 1, 10, $sta_name,$PS_file);
}

# ---- plot T components next -----
print "\n plot T component first, please wait ...\n";

$BSCALE = "-BwSne";
$xshift = $XLEN*1.12;

`gmt basemap -JX$SCALE -R$RSCALE -X$xshift  $BSCALE -Bxa$xanot -Bya$yanot`;
`gmt plot time.of -JX$SCALE -R$RSCALE  -N -W90/90/90 `;
`gmt plot $phase2.of -JX$SCALE -R$RSCALE  -N -W100/100/100`;
&PSTEXT($XMIN +($XMAX-$XMIN)/20,$YMIN+($YMAX-$YMIN)/50,9, 0, 3, 5, "Tangential",$PS_file);
$x = $ph2_time[$YMIN_int]; $y = $YMIN + $y_offset/8;
&PSTEXT($x,$y ,8, $textangle, 3, 5, $phase2,$PS_file);

for($i=0; $i<@T_FILES; $i++)
{
    $Line = "0/0/0";
    $SACFILE = $T_FILES[$i];chomp($SACFILE);
    @tsegs = split(/_/,$SACFILE); 
    $sta_name = substr($tsegs[1],0,5);
    $B_time = $a_times_n[$i] - $X_O;
    $E_time = $B_time + ($XMAX -$XMIN);
    @name = split(/\./,$SACFILE);
    my($SACFILE,$long) = split/\s+/,`saclst kstnm f $SACFILE`; 
    @sta_name = split(/\./,$long);
    #system("$WF $SACFILE $B_time $E_time $XMIN $YMAX $XMAX $YMIN 1 $temp ");
    $scale = $unit_r[$i];
    #print "No.$i, $scale1\n";
    # waveform begin_time end_time infile l_b.x r_t.x height filename
    system("$WFT $SACFILE $B_time $E_time $XMIN $XMAX $scale $temp ");

    `gmt plot $temp -JX$SCALE -R$RSCALE  -N -W$Line`;
    &PSTEXT($XMAX + $x_offset, $gcarc[$i],8, $textangle, 1, 10, $sta_name[0],$PS_file);
}

`gmt basemap -JX$SCALE -R$RSCALE $BSCALE`;
`gmt end`;
`rm $temp`;
chdir $DIR;
#`rm *.BHR *.BHT`;
chdir $Old_dir;


##-----------------------------------------------------------
## subroutines following

sub caltimetable{
    my($phases,$sdep) = @_;
    @phase = split(",",$phases);
    open(TCURVE, "|tcurve");
	for($i = 0; $i < @phase; $i++){
 	    print TCURVE "$phase[$i]\n";
        }
	print TCURVE "\n";
	print TCURVE "$sdep\n";
    close(TCURVE);
}

sub timecurves{
    my($delay,$reduce_vel) = @_;
    open(PT, "< tcurve.out");
	$l = <PT>; chomp($l);
	@phases = split(" ",$l);
	for($i=0; $i < @phases; $i++){
	    $tcurves{$phases[$i]} = "";
        }
	while($l = <PT>){
	    chomp($l);
	    @ts = split(" ",$l);
	    for($i=0; $i < @phases; $i++){
		if($i > 0){
		    $ts[$i] -= $reduce_vel * $ts[0] -$delay;
                }
                # next step using hash: Pn related to arrival times: 22.,33.,46., ...
		$tcurves{$phases[$i]} = "$tcurves{$phases[$i]},$ts[$i]";
            }
        }
    close(PT);
}

sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $textfont, $just, $text, $ps) = @_;
    open(GMT,"| gmt text -N -F+f+a+j ");
        print GMT "$xx $yy $textsize,$textfont $textangle $just $text\n";
    close(GMT);
}
