#####################################################
# perl script for merging the ray model and the prem model into a file
# written by Liang Zhao, IGG CAS
# Feb 18, 2006
#####################################################

$ray_zl = "./ray_zl";

$ray_par = "s51_epi90_grt.par";
$prem_par = "prem.par";
$merge_file = "prem_".$ray_par;

# 1. run the program creating ray model
`$ray_zl`;

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

