#!/bin/bash
#
# Lotka-Volterra demo using Lotka.py:
# compares stochastic with deterministic simulations;
#
# in case of stochastic, due to the small number of molecules used,
# there is a high probability of predator or prey extinction; hence in
# order to produce some nice plots we have to repeat the simulation
# several times, until the ecosystem survives at least for a specified
# minimum amount of time.
#
# by Lidia Yamamoto, Belgium, August 2013
#

python Lotka.py > lotka-det.gr

RETRY=100
MINTIME=30
# retry stoch sim until system survives at least until t=MINTIME
tf=0
i=0
while [ $tf -lt $MINTIME ] ; do
  python Lotka.py stoch > lotka-stoch.gr
  tf=`tail -1 lotka-stoch.gr | cut -f 1 | cut -f 1 -d '.'`
  echo "tf=$tf"
  i=$[$i + 1]
  if [ $i -ge $RETRY ] ; then
    tf=$MINTIME
  fi
done

gnuplot <<EOF
set terminal postscript eps color 20
set ylabel 'Concentration'
set xlabel 'time (s)'
set output 'lotka-det.eps'
plot \
  'lotka-det.gr' using 1:2 w l title 'grass' lt 2 lw 2, \
  'lotka-det.gr' using 1:3 w l title 'rabbit' lt 3 lw 2, \
  'lotka-det.gr' using 1:4 w l title 'fox' lt 1 lw 2 

set output 'lotka-stoch.eps'
plot \
  'lotka-stoch.gr' using 1:2 w l title 'grass' lt 2 lw 2, \
  'lotka-stoch.gr' using 1:3 w l title 'rabbit' lt 3 lw 2, \
  'lotka-stoch.gr' using 1:4 w l title 'fox' lt 1 lw 2
EOF


