#!/opt/antelope/python2.7.2/bin/python2.7
#-----------------------------------------------------
# This script makes a summary of misfit data generated by cap when
# running with only_first_motion=1
#
# Output file can be used by GMT scripts to plot misfit on the lune.
# 
#  usage: % python misfit2lune.py out.misfit.fmp_ > out.misfit.fmp_summary
#
#  input file: out.misfit.fmp_
#       usually large (760MB for 5-deg grid search)
#       example input data: clvd/iso/str/dip/rak/misfit
#
#       -30.0 -71.3   0.0  27.3 -90.0   7
#       -30.0 -71.3  10.0  27.3 -90.0   7
#       -30.0 -71.3  20.0  27.3 -90.0   7
#       -30.0 -71.3  30.0  27.3 -90.0   7
#       ...
#
#   output: summary file out.misfit.fmp_summary
#       gamma
#       delta
#       minimum polarity misfit
#       count of minimum misfit
#       fraction = (misfit0 / misfit_all)
#
#  20130822 cralvizuri
#-----------------------------------------------------

import sys

if len(sys.argv) < 2:
    sys.exit('Usage: %s out.misfit.fmp_ > outputfile ' % sys.argv[0])

inputfile = str(sys.argv[1])
sys.stderr.write('input file: %s \n' % inputfile)

file = open(inputfile, 'r')
#coordinates
c = {}

sys.stderr.write('processing polarities...  ')
# process data summary
for line in file:
    #items
    i = [float(n) for n in line.split()]

    x, y, misfit = i[0], i[1], i[5]

    if (x, y) not in c:
        c[(x, y)] = [misfit, 1, 1]
        ### debug
        #print("{:5.1f} {:5.1f} {} {}".format(
        #    x,y,misfit, c[(x,y)][2]))
    else:   # (x,y) already in c
        c[(x, y)][2] += 1
        if c[(x, y)][0] > misfit:
            # lower misfit found. update
            c[(x, y)][0], c[(x, y)][1] = misfit, 1
        elif c[(x, y)][0] == misfit:
            c[(x, y)][1] += 1

# output results
for key in c:
    print("{:5.1f} {:5.1f} {:4.1f} {} {:.2e}".format(
        key[0], key[1], c[key][0], c[key][1],
        float(c[key][1])/float(c[key][2])))


sys.stderr.write('done.\n')
