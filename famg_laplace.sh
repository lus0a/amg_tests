#!/bin/bash

ugsubmit 1 --- ugshell -ex amg/famg_laplace.lua -numPreRefs 1 -numRefs 4 -amgMaxNodesForBase 2000 -amgXC -outdir $HOME/results
ugsubmit 4 --- ugshell -ex amg/famg_laplace.lua -numPreRefs 1 -numRefs 5 -amgMaxNodesForBase 2000 -amgXC -outdir $HOME/results
ugsubmit 16 --- ugshell -ex amg/famg_laplace.lua -numPreRefs 1 -numRefs 6 -amgMaxNodesForBase 2000 -amgXC -outdir $HOME/results
ugsubmit 64 --- ugshell -ex amg/famg_laplace.lua -numPreRefs 1 -numRefs 7 -amgMaxNodesForBase 2000 -amgXC -outdir $HOME/results
ugsubmit 256 --- ugshell -ex amg/famg_laplace.lua -numPreRefs 1 -numRefs 8 -amgMaxNodesForBase 2000 -amgXC -outdir $HOME/results
