import numpy as np
import math
import os
import sys
import argparse
import math

parser = argparse.ArgumentParser("Do median filtered baseline removal")

parser.add_argument("--decln", type=float, help="Declination", default=0.0)
parser.add_argument("--beamwidth", type=float, help="Beam width", default=1.0)
parser.add_argument("--interval", type=float, help="Sample interval (seconds)", default=20.0)

args = parser.parse_args()

#
# Translate beam-width into transit time at the given declination
#
beamdwidth = args.beamwidth
beamwidth = 4.0 / math.cos(math.radians(args.decln))
beamwidth *= 60.0

# 
# Translate this into the size of the median filter
#
# Rule-of-thumb is this should be 10 times beamwidth
#
MLEN = beamwidth / args.interval
MLEN *= 10
MLEN = int(MLEN)

MLEN=300

marray=[]
mndx = 0

alpha = 1.0/float(MLEN)
beta = 1.0 - alpha
mfilt = None
while True:
    ln = sys.stdin.readline()
    if (len(ln) == 0):
        break
    toks=ln.replace("\n", "")
    toks=ln.split(" ")
    value=float(toks[1])
    if (mndx < MLEN):
        marray.append(value)
        mndx += 1
    else:
        mfilt = float(np.median(marray))
        print ("%.5e" % (value - mfilt))
        marray = [value] + marray[:-1]
