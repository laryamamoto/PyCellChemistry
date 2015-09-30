#!/bin/bash

PREF='matchem'

python MatrixChem.py > ${PREF}.gr

header=`head -1 ${PREF}.gr | cut -f 2-`
pltline='plot '
sep=''
i=1
for clm in $header ; do
  i=$[$i + 1]
  pltline="$pltline$sep\"${PREF}.gr\" using 1:$i w l lw 2 title \"$clm\""
  sep=', '
done

gnuplot <<EOF
set terminal postscript eps color 16
set ylabel 'Molecule count'
set xlabel 'time (s)'
set output "${PREF}.eps"
$pltline
EOF
