#!/bin/bash

PREF='repress'
RESDIR='repress-results'

python Repressilator.py > ${PREF}.gr
head -1 ${PREF}.gr | tr '\t' '\n' > ${PREF}-head.txt

for i in p g c m ; do # p=protein, g=gene, c=complex, m=mRNA
  cmd='plot'
  sep=''
  for j in a b c ; do # 3 variants of each kind, corresponding to the 3 genes
    mol=$i$j
    clm=`grep -n $mol ${PREF}-head.txt | cut -f 1 -d ':'`
    #echo "$i$j = $clm"
    cmd="$cmd$sep \"${PREF}.gr\" using 1:$clm w l lw 2 title \"$mol\""
    sep=','
  done
  #echo "cmd=$cmd"

gnuplot <<EOF
  set terminal postscript eps color 16
  set ylabel 'number of molecules'
  set xlabel 'time (s)'
  set output "${PREF}-${i}.eps"
  $cmd
EOF
done

mkdir $RESDIR
mv ${PREF}*.gr ${PREF}*.eps $RESDIR
\rm -f ${PREF}-head.txt
