# ####################################################
#  perl scripts for filter ncisp-6 sac files, deconvolute 
#  instrument response, 
#  Liang Zhao, IGG CAS
#  Last run 100910
# ####################################################

$DIR    = "/home/zhl/work/s-sks/ncisp6/data/raw";

$RESPONSE = "/home/zhl/work/share/cmg3esp60.response";

chdir $DIR;

@DIRS = `ls `;

for($j = 0; $j < @DIRS; $j++){
   chdir $DIR;
   $dir = $DIRS[$j]; chomp($dir);
   if( -d $dir)   {
      chdir $dir;

#      if(! -d $DIRTMP/$dir) { mkdir("$DIRTMP/$dir", 0744);}
      print stdout "changing to directory $dir\n";

      @STATION = `ls *.sac`;
      for($st = 0; $st < @STATION; $st++ )
      {
         $SACFILE = $STATION[$st]; chomp($SACFILE);
         @segs=split(/\./,$SACFILE);
         @names=split(/_/,$segs[0]);
         $segs[0]=$names[1];
         $newname= join(".",@segs);
         `mv $SACFILE $newname`;
         print STDOUT "the sacfile is $SACFILE\n";	

         open(SAC, "|sac");
         print SAC "r $newname\n";
         print SAC "w over\n";
         print SAC "rmean\n";
         print SAC "rtrend\n";
         print SAC "taper\n";
#         print SAC "bp bu co 0.02 0.5 n 2 p 2\n";    
         print SAC "trans from polezero subtype $RESPONSE to none\n";
         print SAC "bp bu co 0.02 0.5 n 2 p 2\n";
         print SAC "w over\n";
         print SAC "quit\n";
         close(SAC);
  
         print STDOUT "$SACFILE was filtered\n";
      }
   }
}

