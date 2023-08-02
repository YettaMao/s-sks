#!/bin/bash

for ani in {0..180..15}; do
	cp "deprf_${ani}" deprf
	#echo "deprf_${ani}"
	perl reverse.pl
	perl prem_epi83.pl
	cp Fig1.ps "Fig_${ani}.pdf"
done
