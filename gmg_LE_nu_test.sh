#!/bin/bash
oder=0
for nu in 0.0 0.3 0.4 0.45 0.49
do
	for numSmooths in 1 2 4 8 16 32 64 128
	do
		for numRefs in 4 5 6 7 8
		do 
			ugsubmit 1 --- ugshell -ex amg/gmg_LE_nu.lua -numRefs $numRefs -LE -nu $nu -numSmooths $numSmooths -outfile $HOME/results/gmg_LE_res.txt -order $order
			((order++))			
		done
	done
done
