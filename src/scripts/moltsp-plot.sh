#!/bin/bash
#
# moltsp-plot.sh: plots Molecular TSP results; invoked from moltsp.sh
#

NCITIES=10
if [ "$1" != "" ] ; then
  NCITIES=$1
fi

#---------- plot resulting tours using gnuplot

p0=`expr 2 \* $NCITIES / 10 `
p1=`expr 2 \* $NCITIES + $p0 `

for INPF in gen*.gr ; do

PREF=`basename $INPF .gr`
OUTF=${PREF}.eps

gnuplot <<EOF

set terminal postscript eps 20

set size 0.75,1
#set pointsize 2
unset xtic
unset ytic
unset key
set style fill transparent solid 0.7 noborder

set output "$OUTF"

plot [-$p0:+$p1][-$p0:+$p1] \
     "$INPF" using 1:2 with lines, \
     "$INPF" using 1:2 with circles

EOF

done

