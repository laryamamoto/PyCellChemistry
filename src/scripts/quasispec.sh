#!/bin/bash
#
# quasispec.sh: one run of Quasispecies.py, generates outputfiles
#   quasispec-mp*.* with the results of the simulations (time evolution
#   plot; histogram of species distributions at the end of the
#   simulation); see book chapter 7 for more details about the
#   interpretation of the results.
#
# usage: quasispec.sh [ mp ]
#
# where mp is the proportion of mutations with respect to the error threshold
# value (default mp=1.0)
#
# by Lidia Yamamoto, Kraainem, Belgium, July 2013
#

QPREF='quasispec'

mp=1.0
if [ "$1" != "" ] ; then
  mp=$1
fi

function plotqs () {
PREF="$1"
gnuplot <<EOF
  set terminal postscript eps color 20
  set ylabel 'ratio'
  set xlabel 'time (s)'
  set output "${PREF}.eps"
  plot \
    "${PREF}.gr" using 1:2 w l lw 2 title 'relative number of species', \
    "${PREF}.gr" using 1:3 w l lw 2 title 'best fitness', \
    "${PREF}.gr" using 1:5 w l lw 2 title 'average fitness', \
    "${PREF}.gr" using 1:4 w l lw 2 title 'worst fitness'
EOF
}

prefmp="${QPREF}-mp${mp}"
histf="${prefmp}-hist.gr"
python Quasispecies.py $mp $histf > ${prefmp}.gr 2> ${prefmp}.log
plotqs $prefmp
python scripts/histogram.py $histf

