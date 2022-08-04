"""
Embedded Python Blocks:

Each time this file is saved, GRC will instantiate the first class it finds
to get ports and parameters of your block. The arguments to __init__  will
be the parameters. All of them are required to have default values!
"""

import numpy as np
from gnuradio import gr
import time
import ephem


# Given longitude(decimal degrees as a float)
#
# Return the current sidereal time as a string with
#  "," separated tokens
#
def cur_sidereal(longitude):
    longstr = "%02d" % int(longitude)
    longstr = longstr + ":"
    longitude = abs(longitude)
    frac = longitude - int(longitude)
    frac *= 60
    mins = int(frac)
    longstr += "%02d" % mins
    longstr += ":00"
    x = ephem.Observer()
    x.date = ephem.now()
    x.long = longstr
    jdate = ephem.julian_date(x)
    tokens=str(x.sidereal_time()).split(":")
    hours=int(tokens[0])
    minutes=int(tokens[1])
    seconds=int(float(tokens[2]))
    sidt = "%02d,%02d,%02d" % (hours, minutes, seconds)
    return (sidt)


class blk(gr.sync_block):  # other base classes are basic_block, decim_block, interp_block
    """Embedded Python Block example - a simple multiply const"""

    def __init__(self, fft_size=2048,f1=1420.4058e6,f2=1420.4058e6,f3=611e6,f4=1000e6,
        decln=0.0,longitude=-76.03,prefix="./", logtime=10.0):  # only default arguments here
        """arguments to this function show up as parameters in GRC"""
        gr.sync_block.__init__(
            self,
            name='Logger for multi vectors (FFT)',   # will show up in GRC
            in_sig=[(np.float32,fft_size)]*4,
            out_sig=None
        )
        # if an attribute with the same name as a parameter is found,
        # a callback is registered (properties work, too).
        self.fft_size = fft_size
        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
        self.f4 = f4
        self.freqs = [self.f1,self.f2,self.f3,self.f4]
        self.logtime = logtime
        self.fftbufs = []
        for x in range(4):
            self.fftbufs.append(np.zeros(self.fft_size,dtype=np.float32))
        self.fftcounters = [0]*4

        self.decln = decln
        self.prefix = prefix
        self.longitude = longitude
        self.then = time.time()

    def work(self, input_items, output_items):
        #
        # For each input channel
        #
        for x in range(len(input_items)):
            v = input_items[x]
            #
            # For each of the vectors being presented in this channel
            #
            for wectors in v:
                #
                # Integrate
                #
                self.fftbufs[x] = np.add(wectors, self.fftbufs[x])
                self.fftcounters[x] += 1
        #
        # Time to log
        #
        if ((time.time() - self.then) > self.logtime):
            self.then = time.time()
            for x in range(len(self.freqs)):
                if (self.fftcounters[x] > 0):
                    self.fftbufs[x] = np.divide(self.fftbufs[x], self.fftcounters[x])
                    self.fftcounters[x] = 0
                    ltp = time.gmtime()
                    ds = "%04d%02d%02d" % (ltp.tm_year, ltp.tm_mon, ltp.tm_mday)
                    fname = self.prefix+"%04d-%d-" % (int(self.freqs[x]/1.e06), x)
                    fname += ds
                    fname += ".csv"
                    fp = open(fname, "a")
                    # Write some stuff
                    #
                    #
                    curt = cur_sidereal(self.longitude)
                    fp.write("%02d,%02d,%02d" % (ltp.tm_hour, ltp.tm_min, ltp.tm_sec))
                    fp.write(",%s," % curt)
                    fp.write("%.2f," % self.decln)
                    fp.write("%.5f," % (self.freqs[x]/1.0e6))
                    for v in self.fftbufs[x]:
                        fp.write("%.5e," % v)
                    fp.write("\n")
                    self.fftbufs[x] = np.zeros(self.fft_size, dtype=np.float32)
                    fp.close()
                    
                
              
        return len(input_items[0])
