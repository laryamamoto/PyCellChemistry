#!/bin/bash

PREF='dimerstoch'

python DimerStoch.py > ${PREF}.gr
gnuplot <<EOF
set terminal postscript eps color 20
set ylabel 'Concentration'
set xlabel 'time (s)'
set output "${PREF}.eps"

plot \
  "${PREF}.gr" using 1:2 w l lw 2 title 'A', \
  "${PREF}.gr" using 1:3 w l lw 2 title 'B', \
  "${PREF}.gr" using 1:4 w l lw 2 title 'C'
EOF
