#!/bin/bash

PREF='nk'
RESDIR='nk-results'

plotcmd='plot'
sep=''
for N in 2 4 8 16 ; do
  echo "N=$N"
  python NKlandscape.py $N > $PREF$N.gr
  plotcmd="$plotcmd$sep '$PREF$N.gr' using 1:2 w linespoints lw 2 title 'N=$N'"
  sep=','
done

echo 'plotting...'
#echo "$plotcmd" > ${PREF}.gnu
gnuplot <<EOF
set terminal postscript eps color 16
set ylabel 'average distance'
set xlabel 'K'
set output "${PREF}.eps"
$plotcmd
EOF

mkdir $RESDIR
mv ${PREF}*.gr ${PREF}.eps $RESDIR
