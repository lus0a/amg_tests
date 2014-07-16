MYOPTIONS=-walltime 00:20:00

ugsubmit 1 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 3 -numRefs 6 -outdir $HOME/results
ugsubmit 4 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 3 -numRefs 7 -outdir $HOME/results
ugsubmit 16 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 3 -numRefs 8 -outdir $HOME/results
ugsubmit 64 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 3 -numRefs 9 -outdir $HOME/results
ugsubmit 256 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 3 -numRefs 10 -outdir $HOME/results

ugsubmit 1 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 1 -numRefs 1 -outdir $HOME/results -3k -amgXC
ugsubmit 4 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 1 -numRefs 2 -outdir $HOME/results -3k -amgXC
ugsubmit 16 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 1 -numRefs 3 -outdir $HOME/results -3k -amgXC
ugsubmit 64 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 1 -numRefs 4 -outdir $HOME/results -3k -amgXC
ugsubmit 256 $MYOPTIONS --- ugshell -ex amg/famg_laplace_filtering/tv_adaptive_mult.lua -numPreRefs 1 -numRefs 5 -outdir $HOME/results -3k -amgXC


for numVecs in 1 2 3
ugsubmit 1 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 1 -amgXC -3k -numVecs $numVecs -outdir $HOME/results
ugsubmit 4 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 2 -amgXC -3k -numVecs $numVecs -outdir $HOME/results
ugsubmit 16 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 3 -amgXC -3k -numVecs $numVecs -outdir $HOME/results
ugsubmit 64 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 4 -amgXC -3k -numVecs $numVecs -outdir $HOME/results
ugsubmit 256 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 5 -amgXC -3k -numVecs $numVecs -outdir $HOME/results

ugsubmit 1 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 5 -amgXC -numVecs $numVecs -outdir $HOME/results
ugsubmit 4 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 6 -amgXC -numVecs $numVecs -outdir $HOME/results
ugsubmit 16 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 7 -amgXC -numVecs $numVecs -outdir $HOME/results
ugsubmit 64 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 8 -amgXC -numVecs $numVecs -outdir $HOME/results
ugsubmit 256 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 9 -amgXC -numVecs $numVecs -outdir $HOME/results

ugsubmit 1 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 2 -amgXC -numVecs $numVecs -outdir $HOME/results -special
ugsubmit 4 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 3 -amgXC -numVecs $numVecs -outdir $HOME/results -special
ugsubmit 16 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 4 -amgXC -numVecs $numVecs -outdir $HOME/results -special
ugsubmit 64 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 5 -amgXC -numVecs $numVecs -outdir $HOME/results -special
ugsubmit 256 $MYOPTIONS --- ugshell -ex amg/famg_LE_filtering/tv_adaptive_mult.lua -numRefs 6 -amgXC -numVecs $numVecs -outdir $HOME/results -special
end