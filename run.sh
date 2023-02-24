#!/bin/bash
## small script to run the analysis
analysis="MainAnalysis"

##OPTION
echo Which option should I run? 
echo Options are:
echo 0 = run all data and MC one after another
echo 11,12,13,14 = run data only 
echo 2,3,4,5 = run MC samples only ...To be improved
read varname
echo Option is $varname
option=$varname

if ((($option == 11))) ; then
        echo 'WARNING!  Running sample with low mu. There are other you need to add using: hadd data.root dataA.root dataB.root dataC.root dataD.root'
fi
if ((  ($option == 12) || ($option == 13) || ($option == 14) || ($option == 0)  )) ; then
        echo 'WARNING!  Running sample with high mu. There are other you need to add using: hadd data.root dataA.root dataB.root dataC.root dataD.root'
fi

#echo Should I use PROOF? \(will make things faster\)
#echo Options are:
#echo 0 = NO
#echo 1 = YES
#read proofvarname
#echo PROOF option is $proofvarname
#parallel=$proofvarname


## execute and run ROOT
echo "starting ROOT"
##
root -l -b << EOF
.L $analysis.C+
$analysis($option)
EOF
##
echo "end of ROOT execution"
cd /Users/danielernani/Desktop/workspace/MyTopAnalysis/output/  << EOF
EOF