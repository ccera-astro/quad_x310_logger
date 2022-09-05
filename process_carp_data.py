#!/usr/bin/python3
#
#
# Process a set of files that have been produced by the quad_x310_receiver flow-graph
#
#
# Each file contains a collection of records:
#
# UTH,UTM,UTS,LMH,LMM,LMS,DEC,FREQ,BW,[2048 FFT power values]
#
# UTH - UTC Hour
# UTM - UTC Minute
# UTS - UTC Second
# LMH - LMST Hour
# LMM - LMST Minute
# LMS - LMST Second
# DEC - Declination
# FREQ - Center frequency
# BW   - Bandwidth
#
#
# We process each of these records to extra a total-power estimate
#  from the input spectrum.  The spectrum is processed to remove
#  narrow spectral artifacts (both static and dynamic) before
#  computing the total power, and converting into an antenna
#  temperature estimate...
#
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser(description="Process CARP antenna data")

parser.add_argument("file", type=str, help="Input files", metavar="file", nargs="+")
parser.add_argument("--step", type=int, help="Time step", default=10)
parser.add_argument("--tsys", type=float, help="Tsys", default=100.0)
parser.add_argument("--tmin", type=float, help="Dataset minimum temperature", default=15.0)
parser.add_argument("--alpha", type=float, help="Smoothing alpha value", default=0.25)
parser.add_argument("--reduce", type=int, help="Reduce bandwidth amount", default=0)
parser.add_argument("--datastart", type=int, help="Data starting column", default=9)
parser.add_argument("--utc", help="Turn on UTC timestamping", action="store_true", default=False)


args = parser.parse_args()


#
# We keep things "binned" in LMST order
#
# With 86400/args.step bins per sidereal day
#
lmstarray = [0]*int(86400/args.step)
lmstcount = [0]*int(86400/args.step)

for f in args.file:
    sys.stderr.write("Processing %s...\n" % f)
    
    fp = open(f, "r")
    #
    # For each line in the input file
    #
    while True:
        l = fp.readline()
        if (l == ""):
            break
        l = l.replace("\n", "")
        toks = l.split(",")
        
        #
        # Remember the header pieces on each record
        #
        htoks = toks[0:args.datastart]
        
        #
        # The rest are the data
        #
        toks = toks[args.datastart:]
        if (len(toks) != 2049):
            raise ValueError("Input contains incorrect number of tokens")
        
        #
        # Trim off the trailing null token
        #
        if (toks[-1] == ""):
            toks = toks[:-1]
        
        if (args.reduce > 0):
            toks = [0]*args.reduce + toks[args.reduce:-args.reduce] + [0]*args.reduce
        
        #
        # Make into numpy array
        #
        a = np.asarray(toks,dtype=float)
        
        #
        # Reshape to allow 4-point median filter
        #
        a = np.reshape(a,(int(2048/4),4))

        #
        # This will produce a new array of 4-point median values
        #
        # The goal here is to remove "spikes" in the spectrum--likely
        #   RFI.
        #
        m = np.median(a,axis=1)
        
        #
        # Total power is the sum of those values
        #
        s = np.sum(m)
        
        #
        # Determine LMST bin
        #
        tstart=3
        if args.utc == True:
            tstart=0
        lmst = float(htoks[tstart])*3600.0
        lmst += float(htoks[tstart+1])*60.0
        lmst += float(htoks[tstart+2])
        
        lmst = int(lmst)
        lmst = int(lmst / args.step)
        
        lmstarray[lmst] += s
        lmstcount[lmst] += 1
    fp.close()

#
# Determine data minimum
#
minv = 99999.99
for i in range(len(lmstarray)):
    if (lmstcount[i] > 0 and lmstarray[i] > 0):
        s  = lmstarray[i]/lmstcount[i]
        if (s < minv):
            minv = s

#print ("minv %.7f" % minv)
#print ("as temp %.2f" % (minv*(args.tsys+args.tmin)-args.tsys))

#
# Produce a smoothed output
#
aval = -9
a = args.alpha
b = 1.0 - a         
for i in range(len(lmstarray)): 
    if (lmstcount[i] > 0 and lmstarray[i] > 0.0):
        lmst = float(i)*float(args.step)
        lmst /= 3600.0
        
        #
        # Compute average
        #
        s = lmstarray[i]/lmstcount[i]
        
        #
        # Normalize to minv
        #
        t = s/minv
        
        #
        # Convert into (antenna) temperature, given TSYS and TMIN
        #
        t *= (args.tsys+args.tmin)
        t -= args.tsys
        
        #
        # Prime the IIR "pump"
        #
        if (aval < 0):
            aval = t
        
        aval = t*a + b*aval

        print ("%.4f %.6e" % (lmst, aval))
    
