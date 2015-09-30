#!/bin/bash

PREF='logistic'

python Logistic.py > ${PREF}.gr
gnuplot <<EOF
set terminal postscript eps color 20
set ylabel 'concentration of X'
set xlabel 'time (s)'
set output "${PREF}.eps"

plot \
  "${PREF}.gr" using 1:2 w l lw 2 title 'ODE', \
  "${PREF}.gr" using 1:3 w l lw 2 title 'SSA'
EOF

