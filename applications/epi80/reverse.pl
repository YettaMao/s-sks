# perl script to transform semrf -> depthrf

$pi = 3.1415926535897932;
$deg = 6371*$pi/180.;
$top = 446.2710;
$x_o = 8510.3400; # should be consistent with xmin in run.par

$semrf = "semmodel";
$deprf = "deprf";
$hshift = 0; # unit in km

$forward = 1; # from semrf - > depthrf

if($forward ==1) {&forward_transform($semrf,$deprf);}
if($forward ==0) {&backward_transform($semrf,$deprf);}

sub forward_transform{
   my($semrf,$deprf) = @_;
   #$layer: the structural layer number 
   open(FF, "$semrf");
   $l = <FF>; chomp($l);
   $l_num = $l;

   open(DEP,">$deprf");
   print DEP "$l\n";
   for($ll=0;$ll<$l_num;$ll++){ 
      $l  = <FF>; chomp($l);
      print DEP "$l\n";
      $an = <FF>; chomp($an);
      $bn = <FF>; chomp($bn);
      undef(@an_segs);
      undef(@bn_segs);

      @an_segs= split(/ +/,$an);
      @bn_segs= split(/ +/,$bn);
      for($npp=0;$npp<@an_segs;$npp++){
          $ax = ($an_segs[$npp]+$x_o)/$deg;
          printf DEP ("%-11.4f ",$ax);
      }  
      print DEP "\n";
      for($npp=0;$npp<@bn_segs;$npp++){
          $bx = $top - $bn_segs[$npp];
          printf DEP ("%-11.4f ",$bx);
      }  
      print DEP "\n";
   } 
   close(FF);
   close(DEP);
}

sub backward_transform{
   my($semrf,$deprf) = @_;
   #$layer: the structural layer number 
   open(FF, "$deprf");
   $l = <FF>; chomp($l);
   $l_num = $l;

   open(DEP,">$semrf");
   print DEP "$l\n";
   for($ll=0;$ll<$l_num;$ll++){ 
      $l  = <FF>; chomp($l);
      $npts = $segs[0];

      print DEP "$npts  ";
      for($kk=1;$kk<15;$kk++){
          printf DEP ("%-9.4f ",$segs[$kk])
      }
      print DEP "\n";
      $an = <FF>; chomp($an);
      $bn = <FF>; chomp($bn);
      undef(@an_segs);
      undef(@bn_segs);

      @an_segs= split(/ +/,$an);
      @bn_segs= split(/ +/,$bn);
      for($npp=0;$npp<@an_segs;$npp++){
          $ax = $an_segs[$npp]*$deg - $x_o + $hshift;
          printf DEP ("%-11.4f ",$ax);
      }  
      print DEP "\n";
      for($npp=0;$npp<@bn_segs;$npp++){
          $bx = $top - $bn_segs[$npp];
          printf DEP ("%-11.4f ",$bx);
      }  
      print DEP "\n";
   } 
   close(FF);
   close(DEP);
}

