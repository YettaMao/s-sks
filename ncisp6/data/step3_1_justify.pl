# ###################################################
# perl script for justify quality of events one by one
# by view their waveforms
# Run 20140926

$home = "/home/zhl/work/s-sks/ncisp6/data";
$Data = "$home/raw";
$Plot = "plot*/*SKS.ps";
$SUM  = "$home/good_event.txt";

`touch $SUM`;

chdir $Data;
@EVENTS = `ls `;

for($ii=0;$ii<@EVENTS;$ii++)
{
    chdir $Data;
    $event = $EVENTS[$ii]; chomp($event);
    print "change to Event: $event\n";
    if(-d $event)
    {
        chdir $event; 
        open(GSV, "|gs $Plot");
#        print GSV "quit\n";
        $is_good = <>; chomp($is_good);
        if($is_good =~/n/|| $is_good =~/N/)
        {
        }
        else{
            open(SUM,">> $SUM");
            print SUM ("$event\n");
            close(SUM);
        }
        close(GSV);
    }
}

