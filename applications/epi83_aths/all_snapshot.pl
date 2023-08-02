# perl script for plotting all the snapshot
$home = "/home/zhl/work/rf/semsyn";
$snapshot = "$home/plot/snapshot.pl";
$snapshot = "$home/plot/snapshot2.pl";

@outfiles = `ls out*`;
$beg_sec= 480;
$starttime = 300;
$dt = 0.05;
$h = 1.2;
$incre = 160;
@subscripts= ("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z");

for($i=0;$i<@outfiles;$i++){
    $ith = $i +1;
    $fith = $ith;
    if($ith<10) {$fith= "0".$ith;};
    $file = "out"."$ith"; 
    $time = $beg_sec + ($starttime + $i*$incre)*$dt;
    $sub = $subscripts[$i];
    $ps = "snapshot_"."$fith".".ps";
    print "No.$ith perl $snapshot $file $sub $time $h -> $ps\n";
    `perl $snapshot $file $sub $time $h`;
    `mv snapshot.ps $ps`; 
}
