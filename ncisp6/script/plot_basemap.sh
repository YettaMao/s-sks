#gmt: /home/maoyt/anaconda3/envs/maoyt/bin/gmt
# /usr/bin/bash
# Mao 20230228
# Plotting the events distribution of NCISP6 for SKS (80-120)

gmt begin events pdf

#--plot global basemap 
#sta_region=115.5/126.5/39/46.5
stlo=121
stla=43
J=E$stlo/$stla/120/10c
gmt set FORMAT_GEO_MAP=+D
gmt coast -J$J -Rg -A10000 -Ggray -Bg

#--plot stations
gmt plot sta_region.txt -W1.5p,darkblue
gmt plot ncisp6.xy -St0.05c -Gblue

#--plot events
gmt plot event_loc.txt -Sa0.25c -Gred

#--plot circles
echo $stlo $stla 160d |gmt plot -SE- -W1q,red

#--plot text
gmt text -D0c/0.3c << EOF
$stlo -38 80\232 
$stlo -77. 120\232
EOF

gmt end show
