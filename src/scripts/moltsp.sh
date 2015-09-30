#!/bin/bash
#
# moltsp.sh: run MolecularTSP.py, move results to a folder, then plot first
# and last generations
#

NCITIES=10
RESDIR='moltsp-results'
if [ "$1" != "" ] ; then
  NCITIES=$1
fi


#---------- run the MolecularTSP python implementation

python MolecularTSP.py $NCITIES > moltsp.log
if [ $? -ne 0 ]; then
    echo "aborted"
    exit -1
fi

#---------- plot resulting tours using gnuplot

./scripts/moltsp-plot.sh $NCITIES

# move all results to a folder
mkdir $RESDIR
mv moltsp.log gen*.gr gen*.eps $RESDIR
