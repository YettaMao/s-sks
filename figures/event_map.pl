## 
$OUTFILE = "event_map.ps";

`gmtset BASEMAP_TYPE fancy`;
`gmtset ANOT_FONT_SIZE 9`;
`gmtset ANOT_OFFSET 0.04`;
`gmtset FRAME_WIDTH 0.03`;

# Robinson projection
`pscoast -Y15C -X2c -R-75/285/-90/90 -JN105/18 -B60g30 -A10000 -G120 -K -P> $OUTFILE`;
`psxy event_289.xy -R -JN -Sa0.5 -W3/0/0/0 -G255/0/0 -K -O >> $OUTFILE`;
`psxy ncisp6.xy -JN -R -K -O -P -Si0.2 -W5/0/0/220 >> $OUTFILE`;
`psxy gc.xy -JN -M -R -K -O -P  -W5/20/20/20 >> $OUTFILE`;


