# perl script for cut segment of sac files
# Last run on 20140926

$DIR = "/home/zhl/work/s-sks/ncisp6/data/raw";
chdir $DIR;

$pwd = `pwd`; chomp($pwd);
$b = 0;
$e = 1600;

@events = `ls `;

for($ii=0; $ii<@events; $ii++){
   chdir $pwd;
   $event = $events[$ii]; chomp($event);
   if( -d $event){
       print "Enter into $event\n";
       chdir $event;
               open(SAC, "|sac");
               print SAC "cut b $b $e\n";
               print SAC "r *.sac\n";
               print SAC "w over\n"; 
               print SAC "quit\n";
               close(SAC); 
                     
   }
}
