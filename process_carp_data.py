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

from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.units import Quantity
import astropy.constants
import numpy as np
import os
import sys
import argparse
import math
import scipy.signal
import scipy.interpolate
import random
import time


# Direction of motion of the Sun. Info from
# http://herschel.esac.esa.int/hcss-doc-15.0/load/hifi_um/html/hifi_ref_frame.html#hum_lsr_frame
psun = SkyCoord(ra = "18:03:50.29", dec = "+30:00:16.8",frame = "icrs",unit = (u.hourangle,u.deg))
vsun = -20.0*u.km/u.s

# vlsr routine
def vlsr(t,loc,psrc,verbose=False):
    """Compute the line of sight radial velocity

    psrc: SkyCoord object or source
    loc: EarthLocation object of observer
    t: Time object
    """
    # Radial velocity correction to solar system barycenter
    vsrc = psrc.radial_velocity_correction(obstime = t, location = loc)

    # Projection of solar velocity towards the source
    vsun_proj = psrc.cartesian.dot(psun.cartesian)*vsun

    if verbose:
        print("Barycentric radial velocity: {0:+8.3f}".format(vsrc.to(u.km/u.s)))
        print("Projected solar velocity:    {0:+8.3f}".format(vsun_proj.to(u.km/u.s)))
    
    return vsun_proj-vsrc

def doppler_frequency(psrc, t, rest_frequency, loc,  verbose=False):
    """
    Compute the Doppler corrected frequency, taking into account the line of sight radial velocity.

    Args:
        psrc (SkyCoord): sky location for correction
        t (float): time for correction
        rest_frequency (Union[Quantity, float]): observed frequency in LSR
        loc (EarthLocation): observer location

    Returns:
        Quantity: Observable frequency
    """
    v1 = vlsr(t, loc, psrc, verbose=verbose)

    beta = v1/astropy.constants.c
    return np.sqrt((1 + beta)/(1 - beta)) * rest_frequency


def crunchit(indata, crunch):
    out = np.zeros(int(len(indata)/crunch), dtype=np.float64)
    for ndx in range(len(out)):
        for x in range(crunch):
            out[ndx] += indata[(ndx*crunch)+x]
    out = np.divide(out, float(crunch))
    return (out)

def plotspec(fp, indata, freq, bw, scale, offset,sep):
    startf = freq-(bw/2.0)
    incr = bw/len(indata)
    for v in indata:
        fp.write ("%.4f%s%.7e\n" % (startf, sep, (v*scale)+offset))
        startf += incr

def dateificate(s,dateify):
    
    if (dateify == True):
        ltp = time.gmtime(time.time())
        dstr = "%04d%02d%02d" % (ltp.tm_year, ltp.tm_mon, ltp.tm_mday)
        return (s+dstr+"-")
    else:
        return(s)
    
    

parser = argparse.ArgumentParser(description="Process CARP antenna data")

parser.add_argument("file", type=str, help="Input files", metavar="file", nargs="+")
parser.add_argument("--step", type=int, help="Time step (TP only)", default=10)
parser.add_argument("--tsys", type=float, help="System equivalent noise temp (K) (TP Only)", default=100.0)
parser.add_argument("--tmin", type=float, help="Dataset minimum temperature (K) (TP Only)", default=15.0)
parser.add_argument("--alpha", type=float, help="Smoothing alpha value", default=0.25)
parser.add_argument("--reduce", type=int, help="Reduce bandwidth amount (TP Only)", default=0)
parser.add_argument("--datastart", type=int, help="Data starting column", default=9)
parser.add_argument("--utc", help="Turn on UTC timestamping (TP Only)", action="store_true", default=False)
parser.add_argument("--raw", help="Do not convert to Tant (TP Only)", action="store_true", default=False)
parser.add_argument("--db", help="Show as dB above min", action="store_true", default=False)
parser.add_argument("--lmst", help="LMST we're interested in", type=float, default=-1.0)
parser.add_argument("--duration", help="Duration (Hours)", type=float, default=0.0)
parser.add_argument("--fftout", help="FFT output files PREFIX", type=str, default="")
parser.add_argument("--tpout", help="Total power output file PREFIX", type=str, default="")
parser.add_argument("--ctxoffset", help="Sidereal offset for context (minutes) (Spectrral only)", type=float, default=0.0)
parser.add_argument("--redshift", help="Compute red-shift relative to this value (MHz) (Specral only)", type=float, default=0.0)
parser.add_argument("--maskcenter", help="Center frequency for masking (MHz) (TP Only)", type=float, default=0.0)
parser.add_argument("--maskwidth", help="Width for masking (MHz) (TP Only)", type=float, default=0.0)
parser.add_argument("--klen", help="Kernel length for final TP filter", type=int, default=1)
parser.add_argument("--crunch", help="Reduce spectral resolution", type=int, default=int(0))
parser.add_argument("--plotoffset", help="Offset for intermediate plot data (Spectral only)", type=float, default=0.0)
parser.add_argument("--csv", help="Produce CSV files", action="store_true", default=False)
parser.add_argument("--dateify", help="Insert date into filenames", action="store_true", default=False)
parser.add_argument("--vlsr", help="Enable VLSR bin rotation", action="store_true", default=False)
parser.add_argument("--poly", help="Apply 7th-order polynomial fit to baseline data in FFT",
    action="store_true", default=False)
parser.add_argument("--latitude", type=float, help="Local geo latitude", default=45.3491)
parser.add_argument("--longitude", type=float, help="Local geo longitude", default=76.0413)
parser.add_argument("--tra", type=float, help="Target RA as fractional hour-angle", default=5.58333)
parser.add_argument("--tdec", type=float, help="Target DEC as fractional degrees", default=-5.38)
parser.add_argument("--tfreq", type=float, help="Target rest frequency in MHz", default=1424.734)

args = parser.parse_args()

slog = open("shiftlog.dat", "w")

#
# Establish location for astropy routines used later on
#
#
geo_loc = EarthLocation.from_geodetic(lat = args.latitude*u.deg, lon = args.longitude*u.deg, height = 100*u.m)

#
# Establish SkyCoord for target source
#
psrc = SkyCoord(ra = args.tra, dec = args.tdec ,frame = "icrs",unit = (u.hourangle,u.deg))
suffix = ".dat"
sep = " "

if (args.csv):
    suffix = ".csv"
    sep = ","

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

#
# Keep track of record count
#
rcnt = 0
recnum = 0

#
# Used if we're doing VLSR correction
#
# Initial value is 0 and is in units of bins
#
#
shift = 0

for f in args.file:
    sys.stderr.write("Processing %s...\n" % f)
    
    #
    # This "knows" what our filenaming convention is
    #
    filedate = os.path.basename(f)[7:]
    filedate = filedate.split(".")
    filedate = filedate[0]
    
    #sys.stderr.write("The dating is %s...\n" % filedate)
    
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
        
        #
        # Important header data in other parts of the code
        #
        freq = float(htoks[7])
        bw = float(htoks[9])
        decl = float(htoks[6])
        
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
        # We're doing VLSR processing
        # Rotate the array by some number of bins, given by
        #  "shift".
        #
        #
        if (args.vlsr == True):
            l = list(a)
            
            #
            # We don't need to recompute the vlsr correction that often,
            #  and it is rather expensive to compute, so do this only
            #  every 20 records or so
            #
            if ((rcnt % 20) == 0):
                #
                # LMST == RA in our system currently
                #
                
                #
                # Establish pointing to source
                #
                # DECL comes from header in data record
                # Since we're a transit instrument, RA comes from LMST in data record
                #
                #
                
                ras = "%02d:%02d:%02d" % (int(htoks[3]), int(htoks[4]), int(htoks[5]))
                decmin = decl * 60.0
                dech = int(decmin/60.0)
                decmin = abs(decmin) - abs(dech*60.0)
                
                decs = "%02d:%02d" % (dech, decmin)
                #print("decs %s" % decs)

                
                #
                # Pick apart filedate
                #  and reformat into ISO format date
                #
                #
                ts = filedate[0:4]
                ts += "-"
                ts += filedate[4:6]
                ts += "-"
                ts += filedate[6:8]
                
                #
                # Append UTC from input record
                #
                ts = "%sT%02d:%02d:%02d" % (ts, int(htoks[0]), int(htoks[1]), int(htoks[2]) )
                
                #sys.stderr.write("TS: %s\n" % ts)
                
                #
                # Then convert into time acceptable to astropy
                #
                t = Time("%s" % ts, scale="utc",format="isot")
                
                
                #v = vlsr(t,geo_loc,psrc,verbose=False)
                
                #
                # Adjust center frequency based on VLSR
                #
                # Center freq comes from data record header
                #
                #
                fprime = doppler_frequency(psrc, t, args.tfreq, geo_loc)
                
                #
                # Determine difference (in MHz)
                #
                shift = fprime-args.tfreq
                
                #
                # Then determine how many bins that is...
                #
                shift /= (bw/len(a))
                
                #
                # convert to integer--since this is in units of bins
                #
                slog.write("%.5f\n" % shift)
                shift = int(shift)
                shift *= -1

                
                #print ("Updating shift: rcnt %d shift %d fdiff %f ts %s" % (rcnt, shift, fprime-freq, ts))

            #
            # Apply the current shift
            #
            l = l[shift:]+l[:shift]
            a = np.array(l)
        
        #
        # Update our record counter
        #
        rcnt += 1
        
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
    fp = open(dateificate(args.tpout,args.dateify), "w")
    
    values = []
    for v in outbuf:
        values.append(v[1])
    values = np.array(values)
    values = scipy.signal.medfilt(values, kernel_size=args.klen)
    for t,v in zip(outbuf,values):
        fp.write("%.3f%s%.7e\n" % (t[0], sep, v))
        
        
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
    
    #
    # Deal with crunchage
    #
    if (args.crunch > 0):
        ctxarray_low = crunchit(ctxarray_low, args.crunch)
    
    #
    # Record the plot data
    #  
    fp = open(dateificate(args.fftout, args.dateify)+"-context_before%s" % suffix, "w")
    ctxarray_low = np.divide(ctxarray_low, ctxcount_low)
    minv = min(ctxarray_low)
    plotspec(fp, ctxarray_low, freq, bw, 1.0/minv, -args.plotoffset, sep)
    fp.close()
    
    #
    # Deal with crunchage
    #
    if (args.crunch  > 0):
        ctxarray_high = crunchit(ctxarray_high, args.crunch)
    
    #
    # Record the plot data
    #  
    fp = open(dateificate(args.fftout, args.dateify)+"-context_after%s" % suffix, "w")
    ctxarray_high = np.divide(ctxarray_high, ctxcount_high)
    minv = min(ctxarray_high)
    plotspec(fp, ctxarray_high, freq, bw, 1.0/minv, args.plotoffset, sep)
    fp.close()

    
    #
    # "Crunch" the fftarray and ctxarray if indicated with the --crunch
    #  command line option.  This improves sensitivity at the expense of
    #  resolution.
    #
    if (args.crunch > 0):
        fftarray = crunchit(fftarray, args.crunch)
        ctxarray = crunchit(ctxarray, args.crunch)
        
    
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
        
        #
        # Now, apply an aggressive median-filter to that, which
        #  will produce a fairly-smooth baseline estimate
        #
        ctxarray = scipy.signal.medfilt(ctxarray,kernel_size=int(len(ctxarray)/9))
        
        #
        # If "poly", then compute a 7th-order polynomial estimate of
        #  the median filter output, and use that as our final baseline
        #  estimate.  In general, the median-filter output *alone* is
        #  better.
        #
        if (args.poly == True):
            polyfit = np.polyfit(np.arange(0,len(ctxarray)),ctxarray,7)
            p = np.poly1d(polyfit)
            newctx = []
            for x in np.arange(len(ctxarray)):
                newctx.append(p(x))
            ctxarray = newctx
            
        
        #
        # Record the plot data
        #
        fp = open(dateificate(args.fftout, args.dateify)+"-context_merged%s" % suffix, "w")
        plotspec(fp, ctxarray, freq, bw, 1.0, 0.0, sep)
        fp.close()
        
        
        #
        # Normalize the "observation" array
        #
        fftarray = np.divide(fftarray, np.min(fftarray))
        
        #
        # Record the plot data
        #
        fp = open(dateificate(args.fftout, args.dateify)+"-observation%s" % suffix, "w")
        plotspec(fp, fftarray, freq, bw, 1.0, 0.0, sep)
        fp.close()
        
        #
        # Subtract out the "context"
        #
        fftarray = np.subtract(fftarray,ctxarray)
        
        fp = open(dateificate(args.fftout, args.dateify)+"-prescale%s" % suffix, "w")
        plotspec(fp, fftarray, freq, bw, 1.0, 0.0, sep)
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
        fp = open(dateificate(args.fftout, args.dateify)+"-baseline%s" % suffix, "w")
        plotspec(fp, smooth, freq, bw, 1.0, 0.0, sep)
        fp.close()
        
        #
        # Do a little (small kernel) median filtering on the non-smooth version to reduce
        #  narrow RFI blips a bit
        #
        fftarray = scipy.signal.medfilt(fftarray, kernel_size=1)
        
        fp = open(dateificate(args.fftout, args.dateify)+"-observation%s" % suffix, "w")
        plotspec(fp, fftarray, freq, bw, 1.0, 0.0, sep)
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
    fp = open(dateificate(args.fftout, args.dateify)+"-final%s" % suffix, "w")
    
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
            fp.write("%.5f%s%.7e\n" % (freq, sep, smoove))
        else:
            fp.write("%.3f%s%.7e\n" % (rs, sep, smoove))
        freq += incr
    
    
