#!/bin/bash

PREF='chameleon'

python Chameleons.py > ${PREF}.gr
gnuplot <<EOF
set terminal postscript eps color 20
set ylabel 'population'
set xlabel 'simulation time'
set output "${PREF}.eps"

plot \
  "${PREF}.gr" using 1:2 w l lw 2 title 'red', \
  "${PREF}.gr" using 1:3 w l lw 2 title 'green', \
  "${PREF}.gr" using 1:4 w l lw 2 title 'blue'
EOF

