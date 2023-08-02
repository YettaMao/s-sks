# perl script for plotting staggered-grid for anisotropic FD
# Written by Liang Zhao, IGG, CAS
# 20140618
# Modified by Mao, for gmt 6.4.0; 202302

#$Taup="/usr/local/TauP-2.1.1/bin/taup_path";//Zhl
$Taup="/opt/TauP-2.5.0/bin/taup_path";
$PS_file = "syn_s-sks_myt";
$Text = "Texts";
$pi = 3.1415926535;
$ed = 6371;

$xshift = 1;
$yshift  = 3;
$yanot   = 50;
$xanot   = 50;

$max = 1.1;
$XMIN = -$max;
$XMAX = $max;
$YMIN = -$max;
$YMAX = $max;
$YLEN = 20;
$XLEN = 20;

$SCALE  = "$XLEN/$YLEN";
$RSCALE = "$XMIN/$XMAX/$YMIN/$YMAX";
$FSCALE = "-F255/255/255";
$BSCALE = "-Bwsne";
$Line   = "6/0/0/0";

$textsize = 15;
$textangle = 0;

#`gmstset MEASURE_UNIT cm`;
#`gmtset FRAME_PEN 5/255/255/255 `;   

`gmt begin $PS_file pdf`;
#  print "Plotting Frame\n";
`gmt basemap -JX$SCALE -R$RSCALE  -X$xshift -Y$yshift`;
 
$diameter1 = 6371 ;
$diameter2 = 6371 - 2891;
&Draw_arc_rect(-80,80,$diameter1,$diameter2,"arc0.xy"); 
`gmt plot arc0.xy -N -W4p,gray40`;

# CMB anormaly
#$diameter3 = 6371 - 2891 + 400; ;
#$diameter4 = 6371 - 2891;
#&Draw_arc_rect(-35,-10,$diameter3,$diameter4,"arc2.xy"); 
#`gmt plot arc2.xy -N -W1p,red -G150`;

# print "Plotting rays\n";
#@evlos=(-75,-70,-65);
@evlos=(-73,-67,-61);
@stlos=(15,21,27);
$phase= "S";
$phase2= "SKS";
for($ee=0;$ee<@evlos;$ee++){
    $evlo = $evlos[$ee];
    $evla = 0;
    $evdp = 144;
    for($ss=0;$ss<@stlos;$ss++){
    $stlo = $stlos[$ss]; 
    $stla = 0;
# /usr/local/TauP-2.1.1/bin/taup_path -ph P -h 51 -evt 0 -45 -sta 0 0
       &One_ray($evlo, $evla, $evdp, $phase, $stlo, $stla);
       `gmt plot path.of -N -W2.5p,lightbrown`;
       &One_ray($evlo, $evla, $evdp, $phase2, $stlo, $stla);
       `gmt plot path.of -N -W2.5p,green4`;
    }
}

# event
print "output event\n";
&transfer_file("event.xy","temp");
`gmt plot temp -N -Sa0.7 -G255/0/0`;
$bottom = 500;
$diameter2 = 6371 - $bottom;
&Draw_arc_rect(12,27,$diameter1,$diameter2, "arc0.xy"); 
`gmt plot arc0.xy -N -W2p,gray40 -G200/200/200`;

# receiver
print "output receiver\n";
&transfer_file("receiver.xy","temp");
`gmt plot temp -N -Skrtriangle/0.6 -Gblue`;

print "Plotting Texts\n";
&PSTEXTFILE($Text, $PS_file);


`gmt end show`;
`gmtset FRAME_PEN DEFAULT `;

### ########################################################
# for Calculating Ray
# Assign the begining and ending points, and search for a ray connecting them

sub Draw_arc
{
   my($left_a, $right_a, $diameter,$file) = @_;
   $step = 1;
   open(OUT,">$file");
   for($angle = $left_a; $angle<=$right_a; $angle=$angle+$step){
        $ang = $angle/180 * $pi;
        $xx = $diameter/$ed * sin($ang);
        $yy = $diameter/$ed * cos($ang);
        print OUT "$xx $yy\n";
   }
   close(OUT);
}

sub Draw_arc_rect
{
   my($left_a, $right_a, $diameter1, $diameter2, $file) = @_;
   $step = 1;
   open(OUT,">$file");
   for($angle = $left_a; $angle<=$right_a; $angle=$angle+$step){
        $ang = $angle/180. * $pi;
        $xx = $diameter1/$ed * sin($ang);
        $yy = $diameter1/$ed * cos($ang);
        
        if($angle == $left_a){
           $xx0 = $xx;
           $yy0 = $yy;
        }
        print OUT "$xx $yy\n";
   }

   for($angle = $right_a; $angle>=$left_a; $angle=$angle-$step){
        $ang = $angle/180. * $pi;
        $xx = $diameter2/$ed * sin($ang);
        $yy = $diameter2/$ed * cos($ang);
        print OUT "$xx $yy\n";
   }
   print OUT "$xx0 $yy0\n";
   close(OUT);
}

sub transfer_file
{
   my($infile, $file) = @_;
   $step = 1;
   open(INN,"$infile");
   open(OUT,">$file");
   while($line=<INN>){
      chomp($line);
      @segs=split(/ +/,$line);
      $angle = $segs[0];
      $diameter = $segs[1];
      &transfer($angle,$diameter);
      print OUT "$x0 $y0\n";
   }
   close(INN);
   close(OUT);
}


sub One_ray
{
   my($evlo, $evla, $evdp, $phase, $stlo, $stla) = @_;
   print "$Taup -ph $phase -h $evdp -evt $evla $evlo -sta $stla $stlo \n";
   `$Taup -ph $phase -h $evdp -evt $evla $evlo -sta $stla $stlo`;
   &read_TauP_path("taup_path.gmt","path.of");
}

sub read_TauP_path{
   my($in,$file)=@_;
# format 
#    0.00    6320.0      0.00    -45.00
#    0.03    6314.8      0.00    -44.97

   open(IN,$in);
   open(OUT,">$file");
   $line= <IN>; chomp($line);
   $num = 0;
   while($line= <IN>){
      chomp($line);
      @segs= split(/ +/,$line);
      $diam = $segs[2];
      $angle= $segs[4];
      if(!($num > 20 && $diam > (6371-$bottom+20))){
         &transfer($angle,$diam);
         print OUT "$x0 $y0\n";
      }
      $num ++ ;
   }
   close(IN);
   close(OUT);
}

sub transfer{
   my($angle,$diam)= @_;
   $ang = $angle/180 * $pi;
   $x0 = $diam/$ed * sin($ang);
   $y0 = $diam/$ed * cos($ang);   
}

sub PSTEXTFILE{
   my($filename,$PS_file) = @_;
   open(FF, $filename);
   $l = <FF>;
   while ($l = <FF>)
   {
       chomp($l);
       @seg = split(/ +/, $l);
       $x = $seg[0];
       $y = $seg[1];
       $textsize = $seg[2];
       $textangle = $seg[3];
       $title = "";
       for($ii = 4; $ii <@seg; $ii++)
       {
          $title = $title." ".$seg[$ii];
       }
       print "$title\n";
       &PSTEXT($x, $y ,$textsize, $textangle, 10, $title);
   }
   close(FF);
}
sub PSTEXT{
    my($xx, $yy,$textsize, $textangle, $just, $text) = @_;
    open(GMT,"| gmt text -F+f+a+j -R$RSCALE -N -JX$SCALE");
        print GMT "$xx $yy $textsize $textangle $just $text\n";
    close(GMT);
}

