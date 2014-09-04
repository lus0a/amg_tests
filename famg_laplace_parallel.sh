#!/bin/bash

ugsubmit 1 --- ugshell -ex amg/famg_laplace.lua -numRefs 7  -numPreRefs 3 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000
ugsubmit 4 --- ugshell -ex amg/famg_laplace.lua -numRefs 8  -numPreRefs 6 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000
ugsubmit 16 --- ugshell -ex amg/famg_laplace.lua -numRefs 9  -numPreRefs 6 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000
ugsubmit 64 --- ugshell -ex amg/famg_laplace.lua -numRefs 10  -numPreRefs 6 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000
ugsubmit 256 --- ugshell -ex amg/famg_laplace.lua -numRefs 11  -numPreRefs 7 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000
ugsubmit 1024 --- ugshell -ex amg/famg_laplace.lua -numRefs 12  -numPreRefs 8 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000
ugsubmit 4096 --- ugshell -ex amg/famg_laplace.lua -numRefs 13  -numPreRefs 8 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000
ugsubmit 16384 --- ugshell -ex amg/famg_laplace.lua -numRefs 14  -numPreRefs 8 -amgTV lin -amgXC -outfile $HOME/results/laplaceStats.txt -amgMaxNodesForBase 1000

  