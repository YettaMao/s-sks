# perl scripts for plotting the files
# Last run on 20140926

$DIR = "/home/maoyt/work/s-sks/ncisp6/data/raw";

# plot N E components
$PlotSKS     = "/home/maoyt/work/aniso/sws/local/PlotPl/SKS_waveform2.pl ";
# plot R T components with (N,E) files
$PlotSKS_RT  = "/home/maoyt/work/s-sks/ncisp6/data/SKS_waveform_RT.pl";
$WF_C        = "/home/maoyt/work/aniso/sws/local/bin/waveform";

chdir $DIR;

@DIRS = `ls -d *`;

for($j = 0; $j < @DIRS; $j ++)
{
    chdir $DIR;
    $dir = $DIRS[$j]; chomp($dir);
    if( -d $dir){
        chdir $dir;
        @name = split(/_/,$dir);
        $PS_file = "SKS";
        $PS_file = $name[0]."_".$name[1]."_".$name[2]."_".$PS_file;
        print stdout "Plotting directory $dir\n";
        $disp = "plot_".$dir;
#        print "disp in $disp $PS_file \n";
        (-d $disp)||mkdir($disp, 0744);
        if( -d $disp)
        {
             print stdout "changing to directory $disp\n";
             chdir $disp;
             `rm core`;
#             `perl $PlotSKS .. $WF_C $PS_file`;
              print "perl $PlotSKS_RT .. $WF_C $PS_file\n";
	      `perl $PlotSKS_RT .. $WF_C $PS_file`;
        }
    }
}

