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
import scipy.interpolate
import random

def crunchit(indata):
    out = np.zeros(int(len(indata)/2), dtype=np.float64)
    for ndx in range(len(out)):
        out[ndx] = indata[ndx*2]
        out[ndx] += indata[(ndx*2)+1]
    out = np.divide(out, 2.0)
    return (out)

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
parser.add_argument("--duration", help="Duration (Hours)", type=float, default=0.0)
parser.add_argument("--fftout", help="FFT output files PREFIX", type=str, default="")
parser.add_argument("--tpout", help="Total power output file", type=str, default="")
parser.add_argument("--ctxoffset", help="Sidereal offset for context (minutes)", type=float, default=0.0)
parser.add_argument("--redshift", help="Compute red-shift relative to this value (MHz)", type=float, default=0.0)
parser.add_argument("--maskcenter", help="Center frequency for masking (MHz)", type=float, default=0.0)
parser.add_argument("--maskwidth", help="Width for masking (MHz)", type=float, default=0.0)
parser.add_argument("--klen", help="Kernel length for final TP filter", type=int, default=1)
parser.add_argument("--crunch", help="Halve spectral resolution", action="store_true", default=False)


args = parser.parse_args()


FFTSIZE = 2048
C = 299792
DAY = 86400

#
# We keep things "binned" in LMST order
#
# With 86400/args.step bins per sidereal day
#
lmstarray = [0]*int(DAY/args.step)
lmstcount = [0]*int(DAY/args.step)

fftarray = np.zeros(FFTSIZE,dtype=np.float64)
fftcount = 0

ctxarray_low = np.zeros(FFTSIZE,dtype=np.float64)
ctxarray_high = np.zeros(FFTSIZE,dtype=np.float64)

ctxcount_low = 0
ctxcount_high = 0

binwidth = -1

#
# Shortcuts for args.lmst/duration
#
almst = args.lmst * 3600.0
adur = args.duration * 3600.0

recnum = 0
for f in args.file:
    sys.stderr.write("Processing %s...\n" % f)
    
    fp = open(f, "r")
    #
    # For each line in the input file
    #
    while True:
        l = fp.readline()
        recnum += 1
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
        if (len(toks) < FFTSIZE):
            raise ValueError("Input contains incorrect number of tokens")
        
        #
        # Trim off the trailing null token if there is one
        #
        if (toks[-1] == ""):
            toks = toks[:-1]
        
        if (args.reduce > 0):
            toks = [0]*args.reduce + toks[args.reduce:-args.reduce] + [0]*args.reduce
        
        freq = float(htoks[7])
        bw = float(htoks[9])
        
        if (args.maskcenter != 0.0):
            if (binwidth < 0):
                binwidth = bw/FFTSIZE
                startf = freq-bw/2.0
                startmask = args.maskcenter-(args.maskwidth/2.0)
                sndx = (startmask-startf)
                sndx = int(sndx/binwidth)
                endx = sndx + int(args.maskwidth/binwidth)
                
        
        #
        # Make into numpy array
        #
        a = np.asarray(toks,dtype=float)
        
        #
        # Because that implies masking...
        # Apply the mask
        #
        if (binwidth >= 0):
            apwr = np.sum(a)/len(a)
            apwr *= 0.80
            for i in range(sndx,endx):
                a[i] = apwr
        
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
        
        #
        # The actual index value
        # (Which is INT(LMST/stepsize)
        #
        if ((recnum % 2 ) == 0):
            lmsi = math.ceil(lmst / args.step)
        else:
            lmsi = math.floor(lmst / args.step)
        if lmsi >= len(lmstarray):
            lmsi = len(lmstarray)-1

        #
        # We only do "time window" processing if both
        #  duration and lmst are specified on the command line
        #
        if (args.duration > 0.0 and args.lmst >= 0.0):
            
            #
            # If this data item's LMST is within the "interesting window"
            #
            if (lmst >= almst-(adur/2.0) and lmst <= almst+(adur/2.0)):
                lmstarray[lmsi] += s
                lmstcount[lmsi] += 1
                fftarray = np.add(fftarray, a)
                fftcount += 1
            
            #
            # If this data item is within the "surrounding context"
            #   LMST
            #
            
            #
            # First the "below"
            #
            lower = almst-(args.ctxoffset*60.0)
            lower -= adur/2.0
            upper = lower+(adur)
            if (args.ctxoffset > 0.0 and lmst >= lower and lmst <= upper):
                ctxarray_low = np.add(ctxarray_low, a)
                ctxcount_low += 1
            
            #
            # Then the "above"
            #
            lower = almst+(args.ctxoffset*60.0)
            lower -= adur/2.0
            upper = lower+(adur)
            
            #
            # But only record it if they specified an offset for context
            #
            if (args.ctxoffset > 0.0 and lmst >= lower and lmst <= upper):
                ctxarray_high = np.add(ctxarray_high, a)
                ctxcount_high += 1
        else:
            lmstarray[lmsi] += s
            lmstcount[lmsi] += 1
            fftarray = np.add(fftarray, a)
            fftcount += 1
    fp.close()
         
if (args.tpout != "" and args.tpout != None):
    outbuf = []
    
    #
    # Process total-power/continuum data
    #
    outvalues = lmstarray
    outcounts = lmstcount
    for ndx in range(len(outvalues)):
        if (outcounts[ndx] > 0):
            outvalues[ndx] /= outcounts[ndx]
    
    #
    # Determine minv
    #
    minv = 999999.99
    for ndx in range(len(outvalues)):
        if (outcounts[ndx] > 0):
            if (outvalues[ndx] < minv):
                minv = outvalues[ndx]
    #
    # Produce a smoothed output
    #
    aval = -9
    a = args.alpha
    b = 1.0 - a         
    for i in range(len(outvalues)): 
        if outcounts[i] > 0 and outvalues[i] > 0.0:
            lmst = float(i)*float(args.step)
            lmst /= 3600.0
            
            s = outvalues[i]
            
            #
            # Normalize to minv, then use this to
            #  produce a temperature estimate
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
            
            #
            # Convert to dB if requested
            #
            if (args.db == True):
                t = math.log10(t/minv)*10.0
            #
            # Prime the IIR "pump"
            #
            if (aval < 0):
                aval = t
            
            #
            # A single-pole IIR filter
            #
            aval = t*a + b*aval

            outbuf.append((lmst, aval))
    fp = open(args.tpout, "w")
    
    values = []
    for v in outbuf:
        values.append(v[1])
    values = np.array(values)
    values = scipy.signal.medfilt(values, kernel_size=args.klen)
    for t,v in zip(outbuf,values):
        fp.write("%.3f %.5e\n" % (t[0], v))
        
        
    fp.close()

#
# Process FFT data
#
if (args.fftout != "" and args.fftout != ""):
    if (fftcount <= 0):
        raise ValueError("No spectral data within specified range")
    
    
    #
    # Then just smash 'em together for further processing
    #
    ctxarray = np.add(ctxarray_low, ctxarray_high)
    
    
    #
    # We initially decompose into "high" and "low" contexts, just so we can log
    #  them seperately.
    #
    
    if (args.crunch == True):
        ctxarray_low = crunchit(ctxarray_low)
        
    fp = open(args.fftout+"-context_low.dat", "w")
    ctxarray_low = np.divide(ctxarray_low, ctxcount_low)
    minv = min(ctxarray_low)
    for v in ctxarray_low:
        fp.write("%.5e\n" % (v/minv))
    fp.close()
    
    if (args.crunch == True):
        ctxarray_high = crunchit(ctxarray_high)
        
    fp = open(args.fftout+"-context_high.dat", "w")
    ctxarray_high = np.divide(ctxarray_high, ctxcount_high)
    minv = min(ctxarray_high)
    for v in ctxarray_high:
        fp.write("%.5e\n" % (v/minv))
    fp.close()

    
    #
    # "Crunch" the fftarray and ctxarray if indicated with the --crunch
    #  command line option.  This improves sensitivity at the expense of
    #  resolution.
    #
    if (args.crunch == True):
        outarray = np.zeros(int(len(fftarray)/2), dtype=np.float64)
    
        for ndx in range(len(outarray)):
            fndx = ndx * 2
            outarray[ndx] += fftarray[fndx]
            outarray[ndx] += fftarray[fndx+1]
        fftarray  = np.divide(outarray, 2.0)
        
        
        ctxout = np.zeros(int(len(ctxarray)/2), dtype=np.float64)
        for ndx in range(len(ctxout)):
            cndx = ndx * 2
            ctxout[ndx] += ctxarray[cndx]
            ctxout[ndx] += ctxarray[cndx+1]
        ctxarray = np.divide(ctxout, 2.0)
        
    
    #
    # Average all the samples we have in fftarray
    #
    fftarray = np.divide(fftarray, fftcount)
    
    #
    # If there was a context specified, process it
    #
    if (ctxcount_low > 0):
        
        #
        # Rememebr how much the average (baseline) power is
        #  because we're going to normalize both the
        #  observation and the context, and then re-scale
        #  afterwards.
        #
        pwravg = sum(fftarray)
        pwravg /= len(fftarray)
        
        #
        # Determine average of "context" and normalize
        #
        ctxarray = np.divide(ctxarray, ctxcount_low+ctxcount_high)
        ctxarray = np.divide(ctxarray, np.min(ctxarray))
        
        fp = open(args.fftout+"-context.dat", "w")
        for v in ctxarray:
            fp.write("%.5e\n" % v)
        fp.close()
        
        
        #
        # Normalize the "observation" array
        #
        fftarray = np.divide(fftarray, np.min(fftarray))
        
        fp = open(args.fftout+"-observation.dat", "w")
        for v in fftarray:
            fp.write("%.5e\n" % v)
        fp.close()
        
        #
        # Subtract out the "context"
        #
        fftarray = np.subtract(fftarray,ctxarray)
        
        fp = open(args.fftout+"-prescale.dat", "w")
        for v in fftarray:
            fp.write("%.5e\n" % v)
        fp.close()
        
        #
        # Re-scale
        #
        fftarray = np.multiply(fftarray, pwravg)
    
    #
    # Not using "context"--use "self baseline" technique
    #
    if (ctxcount_low <= 0):
        #
        # Compute a large-kernel median filter on observation array, to give a good
        #  baseline estimate
        #
        # It has to be fairly large to "smooth over" RFI blips and spectral
        #   humps.
        #
        smooth = scipy.signal.medfilt(fftarray,kernel_size=177)
        
        #
        # Produce a "baseline" file that shows our self-extracted baseline
        #
        # Useful to help understand the end result
        #
        fp = open(args.fftout+"-baseline.dat", "w")
        for v in smooth:
            fp.write("%.5e\n" % v)
        fp.close()
        
        #
        # Do a little (small kernel) median filtering on the non-smooth version to reduce
        #  narrow RFI blips a bit
        #
        fftarray = scipy.signal.medfilt(fftarray, kernel_size=1)
        
        fp = open(args.fftout+"-observation.dat", "w")
        for v in fftarray:
            fp.write("%.5e\n" % v)
        fp.close()
        
        #
        # Subtract-out the smooth version
        #
        fftarray = np.subtract(fftarray, smooth)
        
        #
        # Add a wee bit to make sure we never take the log of <= 0
        #
        if (args.db == True):
            fftarray = np.add(fftarray,1.0e-5)
            
    #
    # Convert to dB scale
    #
    if (args.db == True):
        fftarray = np.log10(fftarray)
        fftarray = np.multiply(fftarray, 10.0)
        
        #
        # Subtract-out the minimum
        #
        fftarray = np.subtract(fftarray, np.min(fftarray))
    
    #
    # Open the output file
    #
    fp = open(args.fftout+"-final.dat", "w")
    
    #
    # Determine Center frequency and bandwidth from headers
    #
    freq = float(htoks[7])
    bw = float(htoks[9])
    
    if args.redshift > 0.0:
        freq += bw/2.0
        incr = -(bw/len(fftarray))
    else:
        freq -= bw/2.0
        incr = bw/len(fftarray)
        
    a = args.alpha
    b = 1.0 - args.alpha
    smoove = -999.000
    
    #
    # Output the smoothed version
    #
    if (args.redshift > 0.0):
        ls = list(fftarray)
        ls.reverse()
        fftarray = np.array(ls)
    
    for v in fftarray:
        if (args.redshift > 0.0):
            rs = (args.redshift-freq)*1.0e6
            rs /= (args.redshift*1.0e6)
            rs *= C
        if (smoove < -90):
            smoove = v
        smoove = a*v + smoove*b
        if (args.redshift <= 0.0):
            fp.write("%.5f %.5e\n" % (freq, smoove))
        else:
            fp.write("%.3f %.5e\n" % (rs, smoove))
        freq += incr
    
    
