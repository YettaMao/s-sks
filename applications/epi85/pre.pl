#####################################################
# perl script for merging the ray model and the prem model into a file
# written by Liang Zhao, IGG CAS
# June 2, 2006
# Last Run on 20140818
#####################################################

$home = "/home/maoyt/work/s-sks/";

$rayhome   = "/home/maoyt/work/s-sks/";
$ray_zl    = "$rayhome/rays/ray_zl";
$ray_sh_zl = "$rayhome/rays/ray_sh_zl";
$ray_ref   = "$rayhome/ray_ref";

$in = "ray_zl.in";
$sh_in = "ray_sh_zl.in";
$initial = "$home/applications/epi85/initial_model.pl";
$an_model= "$home/applications/epi85/iso2sem.pl";

$ray_par    = "s152_epi85.par";
$ray_sh_par = "s152_epi85_sh.par";
#$prem_par = "/home/zhl/work/aniso/code/rays/prem.par";
$prem_par = "prem_152.par";
$merge_file = "prem_".$ray_par;
$merge_sh_file = "prem_".$ray_sh_par;

$ref = 1; # Using depth reflection phase
$ref = 0; # No using depth reflection phase

# 1. run the program creating ray model
print "1: creating ray model: $ray_par\n";
if(-f $in){
  if($ref) {`$ray_ref`;} 
  else { 
     `$ray_zl`;
  }
}
else{
   print "Error: no P-SV ray parameters read in\n";
}

if(-f $sh_in){
   `$ray_sh_zl`;
}
else{
   print "Error: no SH ray parameters read in\n";
}

# 2. get the source layer number
open(RR,$ray_par);
$line = <RR>;
$line = <RR>; chomp($line);
@seg = split(/ +/,$line);
$source = $seg[1] - 1;
if($seg[0] <100) { $source = $seg[2] -1;}
close(RR);

# 3. read and write prem file for P-SV
print "2. Merge into $merge_file\n";
open(MM, ">$merge_file");
open(PP,"$prem_par");
$line = <PP>; chomp($line);
print MM "$line $source\n";

while($line = <PP>){
   chomp($line);
   print MM "$line\n";
}
close(PP);

open(RR,"$ray_par");
while($line = <RR>){
   chomp($line);
   print MM "$line\n";
}
close(RR);
close(MM); 

# 4. read and write prem file for SH
print "3. Merge into $merge_sh_file\n";
open(MM, ">$merge_sh_file");
open(PP,"$prem_par");
$line = <PP>; chomp($line);
print MM "$line $source\n";

while($line = <PP>){
   chomp($line);
   print MM "$line\n";
}
close(PP);

open(RR,"$ray_sh_par");
while($line = <RR>){
   chomp($line);
   print MM "$line\n";
}
close(RR);
close(MM); 

`rm $ray_par`;
`rm $ray_sh_par`;

# 5. SEM model
print "3: creating SEM model\n";
print "perl $initial\n";
print "perl $an_model\n";
`perl $initial`;
`perl $an_model`;
