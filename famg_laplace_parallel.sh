#!/bin/bash
# parameter 1: additional path/prefix
UGSUBMITARGS=""

#PEMax=`ugsubmit-getconfig pemax`

ugsubmit 1 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 7  -numPreRefs 3 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
ugsubmit 4 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 8  -numPreRefs 6 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
ugsubmit 16 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 9  -numPreRefs 6 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
ugsubmit 64 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 10  -numPreRefs 6 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
ugsubmit 256 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 11  -numPreRefs 7 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000

if [ $UGSUBMIT_TYPE == "Hermit" ] || [ $UGSUBMIT_TYPE == "Hornet" ]; then
	ugsubmit 1024 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 12  -numPreRefs 8 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
	ugsubmit 4096 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 13  -numPreRefs 8 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
	ugsubmit 16384 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 14  -numPreRefs 9 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
	ugsubmit 65536 $UGSUBMITARGS --- ugshell -ex amg/famg_laplace.lua -numRefs 15  -numPreRefs 9 -amgTV lin -amgXC -outfile $HOME/results/$1laplaceStats.txt -amgMaxNodesForBase 4000
fi
  