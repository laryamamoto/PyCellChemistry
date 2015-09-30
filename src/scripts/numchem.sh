#!/bin/bash

PREF='numchem'

#python NumberChem.py > ${PREF}.gr
python NumberChemHO.py > ${PREF}.gr

gnuplot <<EOF
set terminal postscript eps color 20
set ylabel 'molecule count'
set xlabel 'time'
set output "${PREF}.eps"

plot \
  "${PREF}.gr" using 1:2 w l lw 2 title 'number of different species', \
  "${PREF}.gr" using 1:3 w l lw 2 title 'cumulative number of effective collisions', \
  "${PREF}.gr" using 1:4 w l lw 2 title 'number of primes'
EOF

