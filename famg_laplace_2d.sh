#!/bin/bash
# parameter 1: additional path/prefix
UGCALL="ugshell -ex amg/famg_laplace.lua -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats2d.txt -amgMaxNodesForBase 4000"
UGSUBMITARGS=" --- $UGCALL"

index=0

# 4 1024
# 5 4225
# 6 16441
# 7 65536

for numProcs in 1 2 4 8 16 32 64 128 256
do
	for numRefs in 5 6 7 8 9 10
	do 
		ugsubmit $numProcs $UGSUBMITARGS -numRefs $numRefs -numPreRefs 5 -index $index
		index=$((index+1))
	done 
done

