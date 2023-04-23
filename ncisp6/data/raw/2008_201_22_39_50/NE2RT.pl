# perl script for converting ZNE system to ZRT system
# Written by Liang Zhao
# April 28, 2009, In 201 McCone Hall, UCB

# File name in directory
# T_196SG.01.2.2001_241_22_47_24.sac
# -> T_196SG.BHR
@sacns = `ls *.01.2.*.sac`;
$now = `pwd`; chomp($now);
for($ii=0;$ii<@sacns;$ii++){
    $sacn = $sacns[$ii]; chomp($sacn);
    @segs = split(/\./,$sacn);
    $segs[2] = "3"; $sace = join(".",@segs);

# change event depth unit from m to km
# add CMPAZ and CMPINC to the sac header for BH[ZNE] channels
# rmean rtrend
          $evdp = `saclst evdp <$sacn`;
          if($evdp > 1000) {$evdp = $evdp/1000;}
          print "$now: $sacn $sace evdp= $evdp\n";
          open(SAC, "|sac");
          print SAC "r $sacn\n";
          print SAC "ch CMPINC 90.\n";
          print SAC "ch CMPAZ  0.\n";
          print SAC "w over\n";
          print SAC "r $sace\n";
          print SAC "ch CMPINC 90.\n";
          print SAC "ch CMPAZ  90.\n";
          print SAC "w over\n";
          print SAC "r $sacn $sace\n";
          print SAC "ch EVDP $evdp\n";
          print SAC "rmean\n";
          print SAC "rtrend\n";
          print SAC "w over\n";
          print SAC "quit\n";
          close(SAC);  

# rotate (z,n,e) system to (z,r,t) system
          print "\n rotate to zrt system\n";
          $new_time = $segs[3];
          $newr = $segs[0]."."."r";
          $newt = $segs[0]."."."t";

          open(SAC, "|sac");
          print SAC "r $sacn $sace\n";
          print SAC "rotate to gcp normal\n";
          print SAC "w junkr junkt\n";

          print SAC "r junkr\n";
          print SAC "ch KCMPNM BHR\n";
          print SAC "w over\n";
          print SAC "r junkt\n";
          print SAC "ch KCMPNM BHT\n";
          print SAC "w over\n";
          print SAC "quit\n";
          close(SAC);  

# copy to database
          print "\n $newr $newt\n";
          `mv junkr $newr`;
          `mv junkt $newt`;     
}
