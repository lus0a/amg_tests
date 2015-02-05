mypath=`pwd`

CALL="ugshell -ex amg/famg_laplace.lua -outfile $mypath/resUG.txt" 
#CALL="echo $CALL"
CALL="ugsubmit 1 --- $CALL"

index=0

function calc()
{
	$CALL $* -amgTV const -index $index
	index=$((index+1))
	$CALL $* -amgTV lin -index $index
	index=$((index+1))
	$CALL $* -amgTV le -index $index
	index=$((index+1))
	$CALL $* -amgXC -amgTV const -index $index
	index=$((index+1))
	$CALL $* -amgXC -amgTV lin -index $index
	index=$((index+1))
	$CALL $* -amgXC -amgTV le -index $index
	index=$((index+1))
	$CALL $* -amgXC -amgTV lin -amgN2Interpolation -index $index
	index=$((index+1))
}


calc -33k -numRefs 0
calc -3k -numRefs 0
calc -3k -numRefs 1
calc -3k -numRefs 2
calc -cembalo -numRefs 0
calc -gitarre -numRefs 3
calc -gitarre -numRefs 4
calc -isolator -numRefs 4
calc -isolator -numRefs 5
calc -numRefs 7 -distortDomain 0.0
calc -numRefs 7 -distortDomain 0.4

# LE
calc -3k -numRefs 0 -LE -amgDelta 0.9
calc -33k -numRefs 0 -LE -amgDelta 0.9
calc -numRefs 7 -distortDomain 0.0 -LE -amgDelta 0.9
calc -numRefs 7 -distortDomain 0.4 -LE -amgDelta 0.9
