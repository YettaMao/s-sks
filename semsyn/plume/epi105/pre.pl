#####################################################
# perl script for merging the ray model and the prem model into a file
# written by Liang Zhao, IGG CAS
# June 2, 2006
# Last Run on May 25, 2011, Nantes
#       test using only downward phase
#####################################################

$home = "/home/zhl/work/rf/semsyn";

$rayhome= "/home/zhl/work/sem";
$ray_zl = "$rayhome/rays/ray_zl";
$ray_ref = "$rayhome/ray_ref";

$in = "ray_zl.in";
$initial = "$home/plume/epi105/initial_model.pl";
$an_model= "$home/plume/epi105/iso2sem.pl";

$ray_par = "s144_epi105.par";
#$prem_par = "/home/zhl/work/aniso/code/rays/prem.par";
$prem_par = "prem_144.par";
$merge_file = "prem_".$ray_par;

$ref = 1; # Using depth reflection phase
$ref = 0; # No using depth reflection phase

# 1. run the program creating ray model
print "1: creating ray model: $ray_par\n";
if(-f $in){
  if($ref) {`$ray_ref`;} 
  else { `$ray_zl`;}
}
else{
   print "Error: no ray parameters read in\n";
}

# 2. get the source layer number
open(RR,$ray_par);
$line = <RR>;
$line = <RR>; chomp($line);
@seg = split(/ +/,$line);
$source = $seg[1] - 1;
if($seg[0] <100) { $source = $seg[2] -1;}
close(RR);

open(MM, ">$merge_file");

# 3. read and write prem file
print "2. Merge into $merge_file\n";
open(PP,"$prem_par");
$line = <PP>; chomp($line);
print MM "$line $source\n";

while($line = <PP>){
   chomp($line);
   print MM "$line\n";
}
close(PP);

# 4. read and write ray file
open(RR,"$ray_par");
while($line = <RR>){
   chomp($line);
   print MM "$line\n";
}
close(RR);
close(MM); 

`rm $ray_par`;

# 5. SEM model
print "3: creating SEM model\n";
print "perl $initial\n";
print "perl $an_model\n";
`perl $initial`;
`perl $an_model`;
