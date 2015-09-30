#---------------------------------------------------------------------------
#
# histogram.py:
#
# histogram of integer values: on the x-axis we have the integer
# value, on the y-axis we have the number of occurrences of this value
#
# usage: histogram.py <file> [ <barstyle> ] [ <legend> ]
# where file is a tab-separated input file containing the x- and y- values
# (first row of input file is the column label)
#
# by Lidia Yamamoto, Kraainem, Belgium, August 2013

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# list of bar colors and hatches
# barstyle =  0    1    2    3     4    5    6    7       8       9
barcolor = ( 'b', 'g', 'r', 'c',  'm', 'y', 'w', '0.75', 'aqua', 'burlywood' )
barhatch = ( 'x', '/', 'o', '\\', '.', 'x', '*', '/'     'O',    '\\')
# don't work?? '-', '+'

def read_data( fnamein ):
    data = np.loadtxt(fnamein, skiprows=1)
    infd = open(fnamein, 'r')
    xlabel = infd.readline().split('\t')[0]
    infd.close()
    return (xlabel, data)

def plot_data( xlabel, data, fnameout, barstyle, legend ):
    xaxis = data[:,0]  # x-axis is first column in input file
    yaxis = data[:,1]  # y-axis is second column in input file
    maxx = xaxis[-1]   # maximum value found on the x-axis (= last point)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xticks(xaxis)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('number of occurrences')

    width = 1  # assume bar width of 1 (tmp, not generic enough!!)
    offset = 1.0 * width / 2.0  # offset to center each bar on its x value
    xaxis = xaxis - offset # center bar around x value
    c = barcolor[barstyle]
    h = barhatch[barstyle]
    p1 = ax.bar(xaxis, yaxis, width, color=c, hatch=h)
    #ax.bar(xaxis, yaxis, width)
    ax.set_xbound(lower=-offset, upper=maxx)
    if legend != '':
        ax.legend( [ p1 ], [ legend ] )
    plt.savefig(fnameout)
    #plt.show()

usage = 'Usage: python histogram.py <file> [ <barstyle> ] [ <legend> ]'

args = sys.argv[1:]
if (len(args) < 1):
    print >> sys.stderr, usage
    exit(-1)

fnamein = args[0]
(base, ext) = os.path.splitext(fnamein)
fnameout = ("%s.eps" % base)
#print "fnameout=", fnameout

barstyle=0
legend=''

if (len(args) > 1):
    try:
        barstyle=int(args[1])
    except ValueError:
        print "invalid barstyle argument:", args[1]
        print >> sys.stderr, usage
        exit(-1)

if (len(args) > 2):
    legend=args[2]
barstyle = barstyle % len(barcolor)

(xlabel, data) = read_data(fnamein)
plot_data(xlabel, data, fnameout, barstyle, legend)

