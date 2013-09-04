# This script makes a summary of misfit data generated by cap
# when running with only_first_motion=1
# output file can be used by GMT scripts to plot misfit on the lune.
# 
#  usage: % python misfit2lune.py > out.misfit
#
#  input: out.misfit_fm.txt  <-- usually large (760MB for 5-deg grid search)
# example input data: clvd/iso/str/dip/rak/misfit
# ...
# -30.0 -71.3   0.0  27.3 -90.0   7
# -30.0 -71.3  10.0  27.3 -90.0   7
# -30.0 -71.3  20.0  27.3 -90.0   7
# -30.0 -71.3  30.0  27.3 -90.0   7
# ...
#
#  20130822 celso
#-----------------------------------------------------

inputfile="out.misfit_fmp.txt"

file = open(inputfile, 'r')
#coordinates
c = {}

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

# output:
#    key[0] = gamma
#    key[1] = delta
# c[key][0] = minimum misfit
# c[key][1] = count of minimum misfit
# fraction  = (misfit0 / misfit_all)

for key in c:
    print("{:5.1f} {:5.1f} {:4.1f} {} {:.2e}".format(
        key[0], key[1], c[key][0], c[key][1], 
        float(c[key][1])/float(c[key][2])))
