
for numRefs in 3 4 5 6
do
	for nu in 0.1 0.2 0.3 0.4 0.45 0.49 0.499
	do
		ugshell -ex amg/famgLE_1.lua -numRefs $numRefs -nu $nu -bWriteStats
	done
done

exit


maxBase=1000
for stretchX in 0.1 0.5 1.0 2 10 20
do
	for numRefs in 3 4 5
	do
		ugshell -ex amg/famgLEcheck.lua -bWriteStats -writeMatrices  -numRefs $numRefs -maxLevels 8 -maxBase $maxBase -stretchX $stretchX  -delta 0.2 -Hinv 0 -maxParents 2 -threePenalty $threePenalty
		for maxParents in 3 4
		do
			for threePenalty in 1.0 1.5 2.0			 
			do
				ugshell -ex amg/famgLEcheck.lua -bWriteStats -writeMatrices -numRefs $numRefs -maxLevels 8 -maxBase $maxBase -stretchX $stretchX  -delta 0.2 -Hinv 0 -maxParents $maxParents -threePenalty $threePenalty
			done
		done
	done
done


	
for threePenalty in 1.0 1.5 2.0
do
	for maxParents in 2 3 4 
	do
		ugshell -ex amg/famgLEcheck.lua -bWriteStats -writeMatrices  -maxLevels 8 -maxBase $maxBase -g 1 -maxParents $maxParents -threePenalty $threePenalty -delta 0.2
	done
done


