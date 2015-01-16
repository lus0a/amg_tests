
for numRefs in 3 4 5 6 7 8 9
do
	ugshell -ex amg/famgLE_structured_unstretched.lua -maxBase 6000 -numRefs $numRefs -bWriteStats -maxLevels 12
done
