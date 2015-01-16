


delta=0.4


for nu in 0.0 0.3 0.4 0.45 0.49
do
for numRefs in 4 5 6 7
do
for bnd in 1 4 
do
		ugshell -ex amg/famg_LE2d_study1.lua -numRefs $numRefs -maxBase 1000 -maxLevels 5 -delta $delta -maxParents 4 -bnd $bnd -tvLin -stretchX 1 -nu $nu -maxParents 3
		ugshell -ex amg/famg_LE2d_study1.lua -numRefs $numRefs -maxBase 1000 -maxLevels 5 -delta $delta -maxParents 4 -bnd $bnd -tvLE -stretchX 1 -nu $nu -maxParents 3
done
done
done

 