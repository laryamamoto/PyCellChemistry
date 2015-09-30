#!/bin/bash

PREF='disperser'

python Disperser.py > ${PREF}.gr
gnuplot <<EOF
set terminal postscript eps color 20
set ylabel 'number of molecules'
set xlabel 'time (s)'
set output "${PREF}.eps"

plot \
  "${PREF}.gr" using 1:2 w l lw 2 title 'X1', \
  "${PREF}.gr" using 1:3 w l lw 2 title 'X2', \
  "${PREF}.gr" using 1:4 w l lw 2 title 'X3', \
  "${PREF}.gr" using 1:5 w l lw 2 title 'X4'
EOF

