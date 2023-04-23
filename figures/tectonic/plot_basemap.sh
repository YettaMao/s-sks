#!/usr/bin/bash
# gmt 6.0 for plot NCC basemap
# maoyt 20230407

echo Start GMT plotting....
ncchome="/home/maoyt/work/map/NCC/nccb"
bline=30p,100/100/100
bdash=30p,255/255/255,-
fline=30p,20/20/20,-


gmt begin NCCmap pdf
	gmt set MAP_FRAME_TYPE plain
	gmt set FONT 140p,black
	gmt coast -R100/131/29/47 -Jb115.5/36/30/50/16c -Ba -W1/1p,black -S40/120/181 -C154/201/219 -A300
	gmt grdimage @earth_relief_01m -I+d
	gmt plot china-geospatial-data-UTF8/CN-border-La.gmt -W10p

	#fill
	gmt plot $ncchome/east-block.xygeo -G200/200/200@20 -W30p,darkgray 
	gmt plot $ncchome/middle-block.xygeo -G150/150/150@20 -W30p,darkgray
	gmt plot $ncchome/west-block.xygeo -G80/80/80@20 -W30p,darkgray
	gmt coast -R100/131/29/47 -Jb115.5/36/30/50/16c -Ba -W1/1p,black -S40/120/181 -C154/201/219 -A300

	#lines
	gmt plot $ncchome/eastb1.xygeo -W$bline
	gmt plot $ncchome/eastb2.xygeo -W$bline
	gmt plot $ncchome/northb.xygeo -W$bline
	gmt plot $ncchome/southb.xygeo -W$bline
	gmt plot $ncchome/westb.xygeo -W$bline
	gmt plot $ncchome/interial-wb.xygeo -W$bdash
	gmt plot $ncchome/interial-eb.xygeo -W$bdash
	gmt plot $ncchome/TLFT1.xygeo -W$fline
	gmt plot $ncchome/TLFT2.xygeo -W$fline
	gmt plot GGL.xy -W40p,40/120/180 
	gmt plot -W40p,black <<EOF
115.6 39
126.5 39
126.5 46.5
115.6 46.5
115.6 39
EOF

	#labels
	gmt text -F+f+a label.txt 
	
	#ncisp6
	gmt plot ncisp6.xy -St100p -G30/30/30



gmt end show
