#!/bin/bash

for numRefs in 1 2 3 4
do 
	ugshell -ex amg/famg_laplace.lua -numRefs $numRefs -testVectorSmooths 0 -amgTV const -writeMatrices -amgXC -amgMaxNodesForBase 1000 -outfile $HOME/results/laplace_serial3d.txt -dim 3
done

exit 

for numRefs in 4 5 6 7 8 9 10
do 
	ugshell -ex amg/famg_laplace.lua -numRefs $numRefs -testVectorSmooths 0 -amgTV const -stretchX 1 -numPreRefs 4 -writeMatrices -amgXC -amgMaxNodesForBase 1000 -outfile $HOME/results/laplace_serial.txt
done