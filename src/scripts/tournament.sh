#!/bin/bash

PREF='tournament'

python Tournament.py > ${PREF}.gr
gnuplot <<EOF
set terminal postscript eps color 16
set ylabel 'ratio'
set xlabel 'time (s)'
set output "${PREF}.eps"

plot \
  "${PREF}.gr" using 1:2 w l lw 2 title 'number of species over population size', \
  "${PREF}.gr" using 1:3 w l lw 2 title 'best fitness', \
  "${PREF}.gr" using 1:4 w l lw 2 title 'worst fitness'
EOF

