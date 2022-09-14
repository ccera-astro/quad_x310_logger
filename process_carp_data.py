#!/usr/bin/python3
#
#
# Process a set of files that have been produced by the quad_x310_receiver flow-graph
#
#
# Each file contains a collection of records:
#
# UTH,UTM,UTS,LMH,LMM,LMS,DEC,FREQ,BW,CAL,[2048 FFT power values]
#
# UTH - UTC Hour
# UTM - UTC Minute
# UTS - UTC Second
# LMH - LMST Hour
# LMM - LMST Minute
# LMS - LMST Second
# DEC - Declination
# FREQ - Center frequency
# CAL - CAL state
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
import math
import scipy.signal

parser = argparse.ArgumentParser(description="Process CARP antenna data")

parser.add_argument("file", type=str, help="Input files", metavar="file", nargs="+")
parser.add_argument("--step", type=int, help="Time step", default=10)
parser.add_argument("--tsys", type=float, help="Tsys", default=100.0)
parser.add_argument("--tmin", type=float, help="Dataset minimum temperature", default=15.0)
parser.add_argument("--alpha", type=float, help="Smoothing alpha value", default=0.25)
parser.add_argument("--reduce", type=int, help="Reduce bandwidth amount", default=0)
parser.add_argument("--datastart", type=int, help="Data starting column", default=9)
parser.add_argument("--utc", help="Turn on UTC timestamping", action="store_true", default=False)
parser.add_argument("--raw", help="Do not convert to Tant", action="store_true", default=False)
parser.add_argument("--db", help="Show as dB above min", action="store_true", default=False)
parser.add_argument("--lmst", help="LMST we're interested in", type=float, default=-1.0)
parser.add_argument("--duration", help="Duration", type=float, default=0.0)
parser.add_argument("--fftout", help="FFT output file", type=str, default="")
parser.add_argument("--tpout", help="Total power output file", type=str, default="")


args = parser.parse_args()


#
# We keep things "binned" in LMST order
#
# With 86400/args.step bins per sidereal day
#
lmstarray = [0]*int(86400/args.step)
lmstcount = [0]*int(86400/args.step)

fftarray = np.zeros(2048,dtype=np.float64)
fftcount = 0

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
        if (len(toks) < 2048):
            raise ValueError("Input contains incorrect number of tokens")
        
        #
        # Trim off the trailing null token if there is one
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
        # Total power is the sum of those values
        #
        s = np.sum(a)
        
        #
        # Determine LMST bin
        #
        tstart=3
        if args.utc == True:
            tstart=0
        lmst = float(htoks[tstart])*3600.0
        lmst += float(htoks[tstart+1])*60.0
        lmst += float(htoks[tstart+2])
        
        lmsi= int(lmst)
        lmsi = int(lmst / args.step)
        
        almst = args.lmst * 3600.0
        adur = args.duration * 3600.0
        if (args.duration > 0.0 and args.lmst >= 0.0):
            if (lmst >= almst-(adur/2.0) and lmst <= almst+(adur/2.0)):
                lmstarray[lmsi] += s
                lmstcount[lmsi] += 1
                fftarray = np.add(fftarray, a)
                fftcount += 1
        else:
            lmstarray[lmsi] += s
            lmstcount[lmsi] += 1
            fftarray = np.add(fftarray, a)
            fftcount += 1
    fp.close()

if (args.tpout != "" and args.tpout != None):
    fp = open(args.tpout, "w")
    #
    # Determine data minimum
    #
    minv = 99999.99
    for i in range(len(lmstarray)):
        if (lmstcount[i] > 0 and lmstarray[i] > 0):
            s  = lmstarray[i]/lmstcount[i]
            if (s < minv):
                minv = s

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
            if (args.raw == False):
                t = s/minv
                
                #
                # Convert into (antenna) temperature, given TSYS and TMIN
                #
                t *= (args.tsys+args.tmin)
                t -= args.tsys
            else:
                t = s
            
            if (args.db == True):
                t = math.log10(t/minv)*10.0
            #
            # Prime the IIR "pump"
            #
            if (aval < 0):
                aval = t
            
            aval = t*a + b*aval

            fp.write("%.4f %.5f\n" % (lmst, aval))
    
    fp.close()

if (args.fftout != "" and args.fftout != ""):
    if (fftcount <= 0):
        raise ValueError("No spectral data within specified range")
    fftarray = np.divide(fftarray, fftcount)
    smooth = scipy.signal.medfilt(fftarray,kernel_size=77)
    fftarray = scipy.signal.medfilt(fftarray, kernel_size=5)
    fftarray = np.subtract(fftarray, smooth)
    fftarray = np.add(fftarray,1.0e-3)
    if (args.db == True):
        fftarray = np.log10(fftarray)
        fftarray = np.multiply(fftarray, 10.0)
    fftarray = np.subtract(fftarray, np.min(fftarray))
    fp = open(args.fftout, "w")
    freq = float(htoks[7])
    bw = float(htoks[9])
    freq -= bw/2.0
    incr = bw/2048.0
    a = args.alpha
    b = 1.0 - args.alpha
    smoove = -999.000
    for v in fftarray:
        if (smoove < -90):
            smoove = v
        smoove = a*v + smoove*b
        fp.write("%.5f %.5e\n" % (freq, smoove))
        freq += incr
    
    
