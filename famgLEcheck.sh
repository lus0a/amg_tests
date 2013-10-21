# stretchX 0.1
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 0.1  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.5 #---> 0.032
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 0.1  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.5 #---> 0.044

# stretchX 0.5
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 0.5  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.5 #---> 0.21
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 0.5  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.5 #---> 0.21

# stretchX 1.0
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.5 #---> 0.24
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 1.0 #---> 0.22
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.1 -Hinv 0 -maxParents 4 -threePenalty 2.0 #---> 0.17
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 1.0  -delta 0.3 -Hinv 0 -maxParents 2 -threePenalty 2.0 #---> 0.6

# stretchX 2
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 3 -maxLevels 3 -stretchX 2  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 #
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 4 -maxLevels 4 -stretchX 2  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 #
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 2  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 # -- 0.4

# stretchX 10
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 3 -maxLevels 3 -stretchX 10  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 #
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 4 -maxLevels 4 -stretchX 10  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 #
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 10  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 # -- 0.4

# stretchX 20
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 3 -maxLevels 3 -stretchX 20  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 #
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 4 -maxLevels 4 -stretchX 20  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 #
ugshell -ex amg/famgLEcheck.lua -bWriteStats -numRefs 5 -maxLevels 5 -stretchX 20  -delta 0.1 -Hinv 0 -maxParents 3 -threePenalty 1.4 # -- 0.4