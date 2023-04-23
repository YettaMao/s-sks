# perl script for calculating velocity along different direction of propagation
# By Liang Zhao, IGG, CAS
# Jan 23, 2007
# Ref: Crampin, S.,1977, A review of the effects of anisotropic layering on the propagation of seismic waves, Geophys. J. R. astr. Soc., 49, 9-27., Page10

# consider plane wave propagating along the xn axis: uj=aj exp[iw(t-pk xk/c)], j,k = 1,2,3
# eigenvalue proble: (Cinjn - lama*E)a = 0 (i,j=1,2,3)
$n=3; 
#$model = "olivine.ela";
#$model = "TI_my.ela";
$model = "TI_higher.ela";
$model = "TI.ela";
$density = 3.31;

$cal = "./cjcbi0";

# Usage: ./cjcbi0 TI.ela  3 3.3

print "$cal $model $n $density\n";
`$cal $model $n $density`;
