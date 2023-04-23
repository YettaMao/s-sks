#perl script for keeping good events and removing bad events
# Written by Liang Zhao
# Last run on 20140926
#

$home = "/home/zhl/work/s-sks/ncisp6/data";
$DIR  = "$home/raw";
$good = "$home/good_event.txt";

open(GG, "$good");

$num = 0;
while($line =<GG>)
{
   chomp($line);
#   print "a good event is $line \n";
   $events[$num] = $line;
   $num ++;
}

chdir $DIR;
@dirs = `ls `;

for($i=0;$i<@dirs;$i++)
{
   $event = $dirs[$i];chomp($event);
   print "the current event is $event \n";
   $is_good = 0;
   for($nn=0;$nn<@events;$nn++)
   {
       if($event =~ /^$events[$nn]$/) 
       {   
            $is_good = 1;
            print "matched good event is $events[$nn] \n";
       }
   }
   if($is_good == 0) 
   {
       print "$event is a bad event directory, removing it...\n";
       `rm -r $event`;
   }
}
close(GG);
