#!/bin/bash

for ani in {0..180..15}; do
	cp "deprf_${ani}" deprf
	#echo "deprf_${ani}"
	perl reverse.pl
	perl prem_epi83.pl
	cp sks_syn.pdf "sys_syn_${ani}.pdf"
done
