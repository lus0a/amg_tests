
for numRefs in 4 5
do
	ugshell -ex amg/famgLE_structured_unstretched.lua -numRefs $numRefs -maxBase 1000 -delta 0.3 -d 0.0 -XC -check -stretchY 1 -XCgriebel -nu 0.3 -stretchX 1 -writeMatrices -maxLevels 6
done

for numRefs in 4 5
do
	ugshell -ex amg/famgLE_structured_unstretched.lua -numRefs $numRefs -maxBase 1000 -tvLE -delta 0.3 -d 0.0 -XC -check -stretchY 1 -XCgriebel -nu 0.3 -stretchX 1 -writeMatrices -maxLevels 6
done


for numRefs in 4 5
do
	ugshell -ex amg/famgLE_structured_unstretched.lua -numRefs $numRefs -maxBase 1000 -tvLin -delta 0.3 -d 0.0 -XC -check -stretchY 1 -XCgriebel -nu 0.3 -stretchX 1 -writeMatrices -maxLevels 6
done



exit

for nu in 0.0 0.3 0.4
do
	mpirun -np 4 ugshell -ex amg/famgLE_structured_unstretched.lua -numRefs 6 -maxBase 1000 -delta 0.3 -d 0.0 -XC -check -stretchY 1 -XCgriebel -nu $nu -stretchX 1 -writeMatrices -maxLevels 6
done
exit

for nu in 0.0 0.3 0.4
do
	mpirun -np 4 ugshell -ex amg/famgLE_structured_unstretched.lua -numRefs 6 -tvLin -maxBase 1000 -delta 0.3 -d 0.0 -XC -check -stretchY 2 -XCgriebel -nu $nu -stretchX 1 -writeMatrices -maxLevels 6
done

for nu in 0.0 0.3 0.4
do
	mpirun -np 4 ugshell -ex amg/famgLE_structured_unstretched.lua -numRefs 6 -tvLin -maxBase 1000 -delta 0.3 -d 0.0 -XC -check -stretchY 10 -XCgriebel -nu $nu -stretchX 1 -writeMatrices -maxLevels 6
done
