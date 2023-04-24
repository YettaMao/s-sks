# perl script for runing SEM programs

$homedir= "/home/maoyt/work/s-sks/";
$pgsdir = "/home/maoyt/work/s-sks/semsyn";

$psv    = "$pgsdir/bin/aserpsv_sem";
$sh     = "$pgsdir/bin/asersh_sem";
$demult = "$pgsdir/bin/demult";
$demult_sh = "$pgsdir/bin/demult_sh";
$sem       = "$pgsdir/bin/sem2d_psvsh";
#$sem      = "$pgsdir/bin/sem2d";

$plotdir = "$homedir/plot";
$Plot_GRT_Directly= "$plotdir/Plot_grt_directly.pl";
$Plot_GRT_sh_Directly= "$plotdir/Plot_grt_sh_directly.pl";
$read_grt = "$plotdir/read_grt";
$plot_grt = "$plotdir/Plot_GRT.pl";
$plot_grt_sh = "$plotdir/Plot_GRT_sh.pl";

$par    = "run.par";

$GRT = 0;
$GRT_SH = 0;
$plot_result = 1;
$evdp = 492;
$phase = "SKS";

# Aniso=1: zl; other: not invoke
$Aniso   = 1;
$bandpass= 1;
$corner1 = 0.02;
$corner2 = 0.5;
$read_uvw = "$plotdir/read_uvw";
$Tosac   = "$plotdir/uvw2sac";
$sac_plot= "$plotdir/sac_plot";
$plot_uvw = "plot_uvw_sks.pl";
$trace   = 52;
$cal_split="/home/maoyt/work/s-sks/applications/epi83/cal_split.pl";

# step 0: GRT calculations
if($GRT){
   $time = `date`; chomp($time);
   print "GRT begin time: $time\n";
   print "$psv par=$par\n";
   `$psv par=$par >grt.txt`;
 
   print "Plotting GRT Directly => grt_b.ps\n";
   `perl $Plot_GRT_Directly`;

   print "$demult par=$par\n";
   `$demult par=$par`;
   if($plot_result == 1){
      `$read_grt par=$par `;
      `perl $plot_grt`;
   }
}
if($GRT_SH){
   $time = `date`; chomp($time);
   print "GRT begin time: $time\n";
   print "$sh par=$par\n";
   `$sh par=$par >grtsh.txt`;

   print "Plotting GRT SH Directly => grt_b_sh.ps\n";
   `perl $Plot_GRT_sh_Directly`;

   print "$demult_sh par=$par\n";
   `$demult_sh par=$par`;
   if($plot_result == 1){
      `$read_grt par=$par `;
      `perl $plot_grt_sh`;
   }
}

# step 2: SEM calculations
$semprocess = "sem.txt";
if($Aniso){
   `rm out* snap*.ps *.c* ` ;
   $time = `date`; chomp($time);
   print "SEM begin time: $time\n";

   print "$sem par=$par >$semprocess\n";
   `$sem par=$par >$semprocess`;

   $time = `date`; chomp($time);
   print "SEM End time: $time\n";
}

# 3.
&tcurve($evdp,$phase);
print "Begin plot\n";
if($plot_result == 1){
    if($bandpass==1){
       print "$Tosac par=$par \n";
       `$Tosac par=$par  > ./sem_p.txt`;

       open(SAC, "|sac");
       print SAC "r *.cz *.cn *.ce\n";
       print SAC "bp bu co $corner1 $corner2 n 2 p 2\n";
       print SAC "w over\n";
       print SAC "quit\n\n";
       close(SAC);

       print "$sac_plot $trace\n";
#       `$sac_plot $trace 0 0 >> ./sem_p.txt`;
       `$sac_plot $trace >> ./sem_p.txt`;
       print "perl $plot_uvw\n >> ./sem_p.txt\n";
       `perl $plot_uvw`;    
    }
    else{
       `$read_uvw par=$par`;
       `perl $plot_uvw`;    
    }
}

`perl plot_sks.pl`;
#4
print "Call cal_split.pl\n";
#`perl $cal_split`;
$remove_pro_files=1;
if($remove_pro_files == 1){
    `rm *.er  *.fa  *.fi  *.log  *.pma  *.pmi  *.sa  *.si `;
}

sub tcurve{
   my($evdep,$phase) = @_;
   open(TCURVE, "|tcurve");
   print TCURVE "$phase\n";
   print TCURVE "\n";
   print TCURVE "$evdep\n";
   close(TCURVE);
}

sub sks_tab{
   my($gcarc) = @_;
   open(TT,"tcurve.out");
   $l = <TT>; chomp($l);
   @seg = split(" ",$l);
   for($ii=0; $ii<@seg; $ii++){
       if($seg[$ii] =~/^SKSac$/) {$col = $ii;};
   }
#   print "sks at $col column\n";
   $num = 0;
   while($l=<TT>){
      chomp($l);
      @seg = split(" ",$l);
      $g[$num] = $seg[0];
      $sks[$num] = $seg[$col];
#      print "$g[$num] $sks[$num]\n";
      $num ++;
   }
   for($ii=0; $ii<@g-1;$ii++){ 
      if($gcarc>=$g[$ii] && $gcarc< $g[$ii+1]){
          $time = $sks[$ii] + ($gcarc-$g[$ii])/($g[$ii+1]-$g[$ii])*($sks[$ii+1] - $sks[$ii]);
      }
   }
   close(TT);
   return $time;
}

