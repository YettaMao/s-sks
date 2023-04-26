# perl scripts for choosing data from NCISP dataset
# Liang Zhao, IGGCAS
# Oct 19, 2009 
# Last run on 20140926

$DIR[0]= "/media/disk/data_disk/Rawdata/ncisp/ncisp6/events";

$DIRNEW= "/home/maoyt/work/s-sks/ncisp6/data/raw";

$gcarcmin = 80.;
$gcarcmax = 120.;

print stdout "Your choice:$headva $min $max\n";

for($ss= 0;$ss<@DIR;$ss++)
{
   chdir $DIR[$ss];
   print stdout "changing to directory $DIR[$ss]\n";
   @events = `ls `;

   for($j = 0; $j < @events; $j++)
   {
      chdir $DIR[$ss];
      $event = $events[$j]; chomp($event);
      # enter into event directory
      if( -d $event){
          chdir $event;  print stdout "changing to directory $event\n";
          @files = `ls *.sac`;
          $SACFILE=$files[0];chomp($SACFILE);
#         print "No.1 sacfile is $SACFILE \n";
          if(!(-f $SACFILE)) {
               print stdout "the $event is empty!\n";
               chdir("..");
               next;
          } 
          $now=`pwd`;   chomp($now);
          print stdout "the current dir is $now\n";
          for($st = 0; $st < @files; $st++ )
          {
               $SACFILE = $files[$st]; chomp($SACFILE);
               $gcarc=`saclst gcarc < $SACFILE`;chomp($gcarc);
               if($st ==0) {print stdout "$SACFILE gcarc= $gcarc\n";}
               if(($gcarc < $gcarcmax)&&($gcarc > $gcarcmin))
               {
                     @seg = split(/\./,$SACFILE);
#                    @time = split(/_/,$seg[3]);
#                    $event = join(":",@time);
                     $event = $seg[3];
                     chdir($DIRNEW); 
                     if(! -d $event) { mkdir($event,0744);}
                     chdir($now);
                     print stdout "gcarc= $gcarc $event/$SACFILE ---> $DIRNEW/$event/$SACFILE\n";
                     `cp $SACFILE $DIRNEW/$event/$SACFILE`;
                } 
          }  
      }
   }
}
